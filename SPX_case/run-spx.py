# SPDX-FileCopyrightText: 2024 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from __future__ import annotations

from pathlib import Path
import argparse
import sys
import math
import re
import statistics
import openmc

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR if (SCRIPT_DIR / "src").is_dir() else SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_ROOT))

from src.serpent_to_openmc import build_openmc_components
from src.serpent_preprocess import load_serpent_lines
from src.serpent_geom_ir import parse_geometry_ir


def _finite_bounding_box(geometry: openmc.Geometry):
    lower_left, upper_right = geometry.bounding_box
    values = [*lower_left, *upper_right]
    if all(math.isfinite(float(value)) for value in values):
        return tuple(lower_left), tuple(upper_right)
    return None


def _all_finite(lower_left, upper_right) -> bool:
    values = [*lower_left, *upper_right]
    return all(math.isfinite(float(value)) for value in values)


def _merge_bounds(bounds):
    if not bounds:
        return None
    lower_left = [min(pair[0][i] for pair in bounds) for i in range(3)]
    upper_right = [max(pair[1][i] for pair in bounds) for i in range(3)]
    return tuple(lower_left), tuple(upper_right)


def _bounds_from_outside_cells(outside_cells: set[openmc.Cell]):
    bounds = []
    for cell in outside_cells:
        if cell.region is None:
            continue
        try:
            lower_left, upper_right = (~cell.region).bounding_box
        except Exception:
            continue
        if _all_finite(lower_left, upper_right):
            bounds.append((tuple(lower_left), tuple(upper_right)))
    return _merge_bounds(bounds)


def _bounds_from_material_cells(geometry: openmc.Geometry):
    bounds = []
    for cell in geometry.get_all_material_cells().values():
        if cell.region is None:
            continue
        try:
            lower_left, upper_right = cell.region.bounding_box
        except Exception:
            continue
        if _all_finite(lower_left, upper_right):
            bounds.append((tuple(lower_left), tuple(upper_right)))
    return _merge_bounds(bounds)


def _resolve_bounds(
    geometry: openmc.Geometry,
    outside_cells: set[openmc.Cell],
    fallback_lower_left=(-300.0, -300.0, -300.0),
    fallback_upper_right=(300.0, 300.0, 300.0),
):
    geom_bounds = geometry.bounding_box
    outside_bounds = _bounds_from_outside_cells(outside_cells)
    material_bounds = _bounds_from_material_cells(geometry)

    lower_left = []
    upper_right = []
    for axis in range(3):
        lo = float(geom_bounds[0][axis])
        hi = float(geom_bounds[1][axis])
        if not (math.isfinite(lo) and math.isfinite(hi)) and outside_bounds is not None:
            lo = float(outside_bounds[0][axis])
            hi = float(outside_bounds[1][axis])
        if not (math.isfinite(lo) and math.isfinite(hi)) and material_bounds is not None:
            lo = float(material_bounds[0][axis])
            hi = float(material_bounds[1][axis])
        if not (math.isfinite(lo) and math.isfinite(hi)):
            lo = float(fallback_lower_left[axis])
            hi = float(fallback_upper_right[axis])
        lower_left.append(lo)
        upper_right.append(hi)

    if any(upper_right[i] <= lower_left[i] for i in range(3)):
        raise RuntimeError("Failed to determine valid geometry extents.")
    return tuple(lower_left), tuple(upper_right)


def _center_from_bounds(lower_left, upper_right):
    return tuple((lo + hi) * 0.5 for lo, hi in zip(lower_left, upper_right))


def _set_outer_vacuum_boundaries(outside_cells: set[openmc.Cell]) -> int:
    changed = 0
    visited = set()
    for cell in outside_cells:
        if cell.region is None:
            continue
        for surface in cell.region.get_surfaces().values():
            key = id(surface)
            if key in visited:
                continue
            visited.add(key)
            if hasattr(surface, "boundary_type") and surface.boundary_type in (
                None,
                "transmission",
            ):
                surface.boundary_type = "vacuum"
                changed += 1
    return changed


def _create_plots(
    model: openmc.Model,
    lower_left,
    upper_right,
    pixels: int,
    zoom: float = 1.05,
    xy_z: float | None = None,
    xz_y: float | None = None,
) -> None:
    origin = _center_from_bounds(lower_left, upper_right)
    width_x = max(float(upper_right[0] - lower_left[0]) * zoom, 1.0)
    width_y = max(float(upper_right[1] - lower_left[1]) * zoom, 1.0)
    width_z = max(float(upper_right[2] - lower_left[2]) * zoom, 1.0)
    xy_origin = (origin[0], origin[1], float(xy_z) if xy_z is not None else origin[2])
    xz_origin = (origin[0], float(xz_y) if xz_y is not None else origin[1], origin[2])

    plot_xy = openmc.Plot()
    plot_xy.basis = "xy"
    plot_xy.origin = xy_origin
    plot_xy.width = (width_x, width_y)
    plot_xy.pixels = (pixels, pixels)
    plot_xy.color_by = "material"
    plot_xy.filename = "spx_xy"

    plot_xz = openmc.Plot()
    plot_xz.basis = "xz"
    plot_xz.origin = xz_origin
    plot_xz.width = (width_x, width_z)
    plot_xz.pixels = (pixels, pixels)
    plot_xz.color_by = "material"
    plot_xz.filename = "spx_xz"

    model.plots = openmc.Plots([plot_xy, plot_xz])


def _infer_xy_z_from_type9_lattices(input_path: Path, fallback: float) -> float:
    try:
        lines = load_serpent_lines(input_path)
        geom_ir = parse_geometry_ir(lines)
    except Exception:
        return fallback

    mids = []
    lattices = geom_ir.get("lattices", {})
    for record in lattices.values():
        if record.get("lat_type_int") != 9:
            continue

        z_values = list(record.get("z_values", []))
        stack_universes = list(record.get("stack_universes", []))
        if len(z_values) < 2:
            continue

        segment_index = None
        # SPX-like convention: second universe in stack (u_02...) is core section.
        for idx, universe_name in enumerate(stack_universes):
            if idx >= len(z_values) - 1:
                break
            if re.search(r"(?:^|_)02", str(universe_name)):
                segment_index = idx
                break

        if segment_index is None:
            segment_index = min(max(len(stack_universes) // 2, 0), len(z_values) - 2)

        z0 = float(z_values[segment_index])
        z1 = float(z_values[segment_index + 1])
        if math.isfinite(z0) and math.isfinite(z1) and z1 > z0:
            mids.append(0.5 * (z0 + z1))

    if not mids:
        return fallback
    return float(statistics.median(mids))


def _pick_fissionable_cells(
    geometry: openmc.Geometry, sample_count: int = 6
) -> list[openmc.Cell]:
    all_cells = list(geometry.get_all_material_cells().values())
    fissionable = []
    for cell in all_cells:
        material = cell.fill
        if isinstance(material, openmc.Material) and _material_is_fissionable(material):
            fissionable.append(cell)

    # Fall back to any material-filled cells if fissionability cannot be inferred.
    if not fissionable:
        fallback = sorted(all_cells, key=lambda cell: cell.id)
        return fallback[:sample_count]

    fissionable.sort(key=lambda cell: cell.id)
    if len(fissionable) <= sample_count:
        return fissionable
    if sample_count <= 1:
        return [fissionable[0]]

    max_index = len(fissionable) - 1
    indices = sorted(
        {
            int(round(i * max_index / (sample_count - 1)))
            for i in range(sample_count)
        }
    )
    return [fissionable[index] for index in indices]


def _material_is_fissionable(material: openmc.Material) -> bool:
    # Avoid hasattr() on properties, because it may execute property code and raise.
    try:
        value = getattr(material, "fissionable")
    except Exception:
        value = None
    else:
        try:
            return bool(value)
        except Exception:
            pass

    fissile_elements = {"U", "Np", "Pu", "Am", "Cm", "Th"}
    entries = getattr(material, "nuclides", None) or getattr(material, "_nuclides", [])
    for entry in entries:
        name = getattr(entry, "name", None)
        if name is None and isinstance(entry, tuple) and entry:
            name = entry[0]
        if name is None:
            continue

        token = str(name).replace("-", "")
        match = re.match(r"^([A-Z][a-z]?)(\d+)", token)
        if not match:
            continue
        element = match.group(1)
        mass_number = int(match.group(2))
        if element in fissile_elements and mass_number >= 230:
            return True

    return False


def _configure_tallies(model: openmc.Model, lower_left, upper_right) -> None:
    tallies = openmc.Tallies()

    reaction_scores = [
        "kappa-fission",
        "fission",
        "nu-fission",
        "absorption",
        "flux",
    ]

    global_power = openmc.Tally(name="global_power")
    global_power.scores = ["kappa-fission", "fission"]
    tallies.append(global_power)

    global_reactions = openmc.Tally(name="global_reaction_rates")
    global_reactions.scores = reaction_scores
    tallies.append(global_reactions)

    kinetics_ifp = openmc.Tally(name="kinetics_ifp")
    kinetics_ifp.scores = [
        "ifp-time-numerator",
        "ifp-beta-numerator",
        "ifp-denominator",
    ]
    tallies.append(kinetics_ifp)

    mesh_xy = openmc.RegularMesh(name="xy_power_mesh")
    mesh_xy.lower_left = lower_left
    mesh_xy.upper_right = upper_right
    mesh_xy.dimension = (20, 20, 1)

    xy_power = openmc.Tally(name="xy_power_map")
    xy_power.filters = [openmc.MeshFilter(mesh_xy)]
    xy_power.scores = ["kappa-fission"]
    tallies.append(xy_power)

    xy_reactions = openmc.Tally(name="xy_reaction_map")
    xy_reactions.filters = [openmc.MeshFilter(mesh_xy)]
    xy_reactions.scores = reaction_scores
    tallies.append(xy_reactions)

    mesh_xz = openmc.RegularMesh(name="xz_power_mesh")
    mesh_xz.lower_left = lower_left
    mesh_xz.upper_right = upper_right
    mesh_xz.dimension = (20, 1, 20)

    xz_power = openmc.Tally(name="xz_power_map")
    xz_power.filters = [openmc.MeshFilter(mesh_xz)]
    xz_power.scores = ["kappa-fission"]
    tallies.append(xz_power)

    xz_reactions = openmc.Tally(name="xz_reaction_map")
    xz_reactions.filters = [openmc.MeshFilter(mesh_xz)]
    xz_reactions.scores = reaction_scores
    tallies.append(xz_reactions)

    for cell in _pick_fissionable_cells(model.geometry, sample_count=6):
        tally = openmc.Tally(name=f"cell_power_{cell.id}")
        tally.filters = [openmc.CellFilter([cell.id])]
        tally.scores = ["kappa-fission", "fission", "nu-fission", "absorption", "flux"]
        tallies.append(tally)

    model.tallies = tallies


def _configure_settings(
    model: openmc.Model,
    lower_left,
    upper_right,
    particles: int,
    batches: int,
    inactive: int,
    temperature_method: str = "interpolation",
    ifp_n_generation: int = 10,
) -> None:
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    if ifp_n_generation > 0:
        settings.ifp_n_generation = min(ifp_n_generation, inactive)
    settings.temperature = {
        "method": temperature_method,
        "range": (250.0, 2500.0),
    }
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
    )
    model.settings = settings


def build_test_model(
    input_path: Path,
    particles: int,
    batches: int,
    inactive: int,
    temperature_method: str = "interpolation",
    ifp_n_generation: int = 10,
    plot_zoom: float = 1.05,
    xy_z: float | None = None,
    xz_y: float | None = None,
) -> openmc.Model:
    materials, geom_components, root = build_openmc_components(input_path)

    boundaries_changed = _set_outer_vacuum_boundaries(geom_components["outside_cells"])
    print("vacuum boundary surfaces set:", boundaries_changed)

    model = openmc.Model(geometry=openmc.Geometry(root))
    model.materials = openmc.Materials(materials.values())

    bbox_direct = _finite_bounding_box(model.geometry)
    if bbox_direct is None:
        print("geometry bounding_box is non-finite; using fallback bounds resolution")
        lower_left, upper_right = _resolve_bounds(
            model.geometry, geom_components["outside_cells"]
        )
    else:
        lower_left, upper_right = bbox_direct

    # Use the finite "inside-core" bounds (from outside cells) for clearer plots.
    plot_bounds = _bounds_from_outside_cells(geom_components["outside_cells"])
    if plot_bounds is None:
        plot_bounds = _bounds_from_material_cells(model.geometry)
    if plot_bounds is None:
        plot_bounds = (lower_left, upper_right)
    plot_lower_left, plot_upper_right = plot_bounds
    default_xy_z = 0.5 * (float(plot_lower_left[2]) + float(plot_upper_right[2]))
    if xy_z is None:
        xy_z = _infer_xy_z_from_type9_lattices(input_path, default_xy_z)

    _create_plots(
        model,
        plot_lower_left,
        plot_upper_right,
        pixels=1000,
        zoom=plot_zoom,
        xy_z=xy_z,
        xz_y=xz_y,
    )
    _configure_tallies(model, lower_left, upper_right)
    _configure_settings(
        model,
        lower_left,
        upper_right,
        particles,
        batches,
        inactive,
        temperature_method=temperature_method,
        ifp_n_generation=ifp_n_generation,
    )

    geometry = model.geometry
    print("materials:", len(model.materials))
    print("cells:", len(geometry.get_all_cells()))
    print("universes:", len(geometry.get_all_universes()))
    print("lattices:", len(geometry.get_all_lattices()))
    print(
        "bbox:",
        tuple(round(float(v), 3) for v in lower_left),
        tuple(round(float(v), 3) for v in upper_right),
    )
    print(
        "plot bbox:",
        tuple(round(float(v), 3) for v in plot_lower_left),
        tuple(round(float(v), 3) for v in plot_upper_right),
    )
    print("plot xy z:", round(float(xy_z), 3))
    return model


def _latest_statepoint(output_dir: Path):
    statepoints = sorted(output_dir.glob("statepoint.*.h5"))
    return statepoints[-1] if statepoints else None

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Parse Serpent SPX to OpenMC, generate xy/xz plots, add power tallies, "
            "and run a low-particle eigenvalue simulation."
        )
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        type=Path,
        default=SCRIPT_DIR / "SPX",
        help="Serpent input file path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory for XML, plots, and statepoint files.",
    )
    parser.add_argument(
        "--particles",
        type=int,
        default=10000,
        help="Particles per generation.",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=60,
        help="Total batches.",
    )
    parser.add_argument(
        "--inactive",
        type=int,
        default=15,
        help="Inactive batches.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="OpenMC threads.",
    )
    parser.add_argument(
        "--temperature-method",
        choices=("interpolation", "nearest"),
        default="interpolation",
        help="Cross-section temperature treatment in OpenMC.",
    )
    parser.add_argument(
        "--ifp-n-generation",
        type=int,
        default=10,
        help=(
            "IFP ancestor generations for kinetics tallies (must be <= inactive). "
            "Set 0 to disable IFP storage."
        ),
    )
    parser.add_argument(
        "--plot-zoom",
        type=float,
        default=1.05,
        help="Multiplier for XY/XZ plot width around the core envelope.",
    )
    parser.add_argument(
        "--xy-z",
        type=float,
        default=None,
        help="Optional z-plane for XY plot (cm).",
    )
    parser.add_argument(
        "--xz-y",
        type=float,
        default=None,
        help="Optional y-plane for XZ plot (cm).",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Skip geometry plotting.",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Skip Monte Carlo run.",
    )
    args = parser.parse_args()

    model = build_test_model(
        args.input_file,
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
        temperature_method=args.temperature_method,
        ifp_n_generation=args.ifp_n_generation,
        plot_zoom=args.plot_zoom,
        xy_z=args.xy_z,
        xz_y=args.xz_y,
    )

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    model.export_to_xml(directory=output_dir)
    model.export_to_model_xml(path=output_dir / "model.xml")

    if not args.skip_plot:
        model.plot_geometry(
            cwd=output_dir,
            directory=output_dir,
            export_model_xml=False,
        )
        print("plots written:", output_dir / "spx_xy.png", output_dir / "spx_xz.png")

    if not args.skip_run:
        openmc.run(cwd=output_dir, threads=args.threads)
        statepoint = _latest_statepoint(output_dir)
        if statepoint is not None:
            with openmc.StatePoint(statepoint) as sp:
                print("statepoint:", statepoint.name)
                print("k-effective:", sp.keff)
        else:
            print("no statepoint file found")


if __name__ == "__main__":
    main()
