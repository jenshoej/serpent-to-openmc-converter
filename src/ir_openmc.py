# SPDX-FileCopyrightText: 2024 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from __future__ import annotations

from math import sqrt
from collections import Counter
from typing import Any, Dict, Iterable, List, Tuple

import openmc
import openmc.data
from openmc.model.surface_composite import CompositeSurface

from .serpent_geometry import sqc, vertical_stack


class HexagonalPrismZ(CompositeSurface):
    """Finite hexagonal prism composed of planar surfaces and z planes."""

    _surface_names = ("prism", "zmin", "zmax")

    def __init__(
        self,
        edge_length: float,
        orientation: str,
        origin: tuple[float, float],
        zmin: float,
        zmax: float,
        boundary_type: str = "transmission",
    ) -> None:
        self.prism = openmc.model.HexagonalPrism(
            edge_length=edge_length,
            orientation=orientation,
            origin=origin,
            boundary_type=boundary_type,
        )
        self.zmin = openmc.ZPlane(z0=zmin, boundary_type=boundary_type)
        self.zmax = openmc.ZPlane(z0=zmax, boundary_type=boundary_type)

    def __neg__(self):
        return -self.prism & +self.zmin & -self.zmax

    def __pos__(self):
        return +self.prism | -self.zmin | +self.zmax


def iter_natural_isotopes(element: str) -> List[tuple[str, float]]:
    abundances = openmc.data.NATURAL_ABUNDANCE.get(element)
    if abundances is None:
        return []

    isotopes: List[tuple[str, float]] = []
    if isinstance(abundances, dict):
        items = abundances.items()
    else:
        items = abundances
    for key, value in items:
        if isinstance(key, int):
            name = f"{element}{key}"
        else:
            name = str(key).replace("-", "")
        isotopes.append((name, float(value)))
    return isotopes


def add_element_expanded(
    material: openmc.Material,
    element: str,
    percent: float,
    mode: str,
    expand: bool,
) -> None:
    if not expand:
        material.add_element(element, percent, mode)
        return

    isotopes = iter_natural_isotopes(element)
    if not isotopes:
        material.add_element(element, percent, mode)
        return

    if mode == "ao":
        for name, frac in isotopes:
            material.add_nuclide(name, percent * frac, mode)
        return

    masses = [openmc.data.atomic_mass(name) for name, _ in isotopes]
    weight_sum = sum(frac * mass for (name, frac), mass in zip(isotopes, masses))
    for (name, frac), mass in zip(isotopes, masses):
        weight_frac = frac * mass / weight_sum
        material.add_nuclide(name, percent * weight_frac, mode)


def materials_from_ir(
    materials_ir: Dict[str, Dict[str, Any]],
    expand_natural_elements: bool = True,
) -> Dict[str, openmc.Material]:
    """Convert material/mix IR records into OpenMC Material objects."""
    openmc_materials: Dict[str, openmc.Material] = {}

    for name, record in materials_ir.items():
        if record["type"] == "material":
            mat = openmc.Material(name=name, material_id=record["material_id"])
            if record["density_mode"] == "mass":
                mat.set_density("g/cm3", record["density"])
            elif record["density_mode"] == "atom":
                mat.set_density("atom/b-cm", record["density"])

            if record["temperature"] is not None:
                mat.temperature = record["temperature"]

            for sab in record["sab"]:
                mat.add_s_alpha_beta(sab)

            for nuclide in record["nuclides"]:
                percent = abs(nuclide["percent"])
                mode = nuclide["percent_mode"]
                if nuclide["A"] is None or nuclide["A"] == 0:
                    add_element_expanded(
                        mat,
                        nuclide["element"],
                        percent,
                        mode,
                        expand=expand_natural_elements,
                    )
                else:
                    mat.add_nuclide(nuclide["name"], percent, mode)

            openmc_materials[name] = mat
        elif record["type"] == "mix":
            mix = [openmc_materials[comp["material"]] for comp in record["components"]]
            mix_per = [comp["percent"] for comp in record["components"]]
            mixed = openmc.Material.mix_materials(mix, mix_per, record["mix_type"])
            if record["material_id"] is not None:
                mixed.id = record["material_id"]
            mixed.name = name
            openmc_materials[name] = mixed

    return openmc_materials


def surfaces_from_ir(
    surfaces_ir: Dict[str, Dict[str, Any]]
) -> Dict[str, openmc.Surface]:
    """Convert surface IR records into OpenMC Surface/CompositeSurface objects."""
    openmc_surfaces: Dict[str, openmc.Surface] = {}

    for name, record in surfaces_ir.items():
        surface_type = record["surface_type"].lower()
        boundary_type = record.get("boundary_type")
        coefficients = list(record.get("coefficients", []))
        uid = int(name) if name.isnumeric() else None
        kwargs: Dict[str, Any] = {"name": name}
        if uid is not None:
            kwargs["surface_id"] = uid
        if boundary_type:
            kwargs["boundary_type"] = boundary_type

        if surface_type == "px":
            surface = openmc.XPlane(coefficients[0], **kwargs)
        elif surface_type == "py":
            surface = openmc.YPlane(coefficients[0], **kwargs)
        elif surface_type == "pz":
            surface = openmc.ZPlane(coefficients[0], **kwargs)
        elif surface_type in ("cyl", "cylz"):
            if len(coefficients) == 3:
                x0, y0, r = coefficients
                surface = openmc.ZCylinder(x0, y0, r, **kwargs)
            elif len(coefficients) == 5:
                x0, y0, r, z0, z1 = coefficients
                center_base = (x0, y0, z0)
                height = z1 - z0
                surface = openmc.model.RightCircularCylinder(
                    center_base, height, r, axis="z"
                )
            else:
                raise ValueError(f"Surface '{name}' has invalid cyl parameters.")
        elif surface_type == "cylx":
            if len(coefficients) == 3:
                y0, z0, r = coefficients
                surface = openmc.XCylinder(y0, z0, r, **kwargs)
            elif len(coefficients) == 5:
                y0, z0, r, x0, x1 = coefficients
                center_base = (x0, y0, z0)
                height = x1 - x0
                surface = openmc.model.RightCircularCylinder(
                    center_base, height, r, axis="x"
                )
            else:
                raise ValueError(f"Surface '{name}' has invalid cylx parameters.")
        elif surface_type == "cyly":
            if len(coefficients) == 3:
                x0, z0, r = coefficients
                surface = openmc.YCylinder(x0, z0, r, **kwargs)
            elif len(coefficients) == 5:
                x0, z0, r, y0, y1 = coefficients
                center_base = (x0, y0, z0)
                height = y1 - y0
                surface = openmc.model.RightCircularCylinder(
                    center_base, height, r, axis="y"
                )
            else:
                raise ValueError(f"Surface '{name}' has invalid cyly parameters.")
        elif surface_type == "sqc":
            x0, y0, half_width = coefficients
            surface = sqc(x0, y0, half_width)
        elif surface_type == "torx":
            x0, y0, z0, a, b, c = coefficients
            surface = openmc.XTorus(x0, y0, z0, a, b, c, **kwargs)
        elif surface_type == "tory":
            x0, y0, z0, a, b, c = coefficients
            surface = openmc.YTorus(x0, y0, z0, a, b, c, **kwargs)
        elif surface_type == "torz":
            x0, y0, z0, a, b, c = coefficients
            surface = openmc.ZTorus(x0, y0, z0, a, b, c, **kwargs)
        elif surface_type == "sph":
            x0, y0, z0, r = coefficients
            surface = openmc.Sphere(x0, y0, z0, r, **kwargs)
        elif surface_type == "plane":
            a, b, c, d = coefficients
            surface = openmc.Plane(a, b, c, d, **kwargs)
        elif surface_type == "cone":
            x0, y0, z0, r, h = coefficients
            center_base = (x0, y0, z0)
            surface = openmc.model.ConicalFrustum(center_base, (0.0, 0.0, h), r, 0.0)
        elif surface_type in ("hexyc", "hexxc"):
            if len(coefficients) < 3:
                raise ValueError(f"Surface '{name}' missing parameters for {surface_type}.")
            x0, y0, d = coefficients[:3]
            edge_length = 2.0 * d / sqrt(3.0)
            orientation = "x" if surface_type == "hexyc" else "y"
            surface = openmc.model.HexagonalPrism(
                edge_length=edge_length,
                orientation=orientation,
                origin=(x0, y0),
                boundary_type=boundary_type or "transmission",
            )
        elif surface_type in ("hexxprism", "hexyprism"):
            if len(coefficients) != 5:
                raise ValueError(
                    f"Surface '{name}' missing parameters for {surface_type}."
                )
            x0, y0, d, zmin, zmax = coefficients
            edge_length = 2.0 * d / sqrt(3.0)
            orientation = "x" if surface_type == "hexyprism" else "y"
            surface = HexagonalPrismZ(
                edge_length=edge_length,
                orientation=orientation,
                origin=(x0, y0),
                zmin=zmin,
                zmax=zmax,
                boundary_type=boundary_type or "transmission",
            )
        elif surface_type == "cuboid":
            xmin, xmax, ymin, ymax, zmin, zmax = coefficients
            surface = openmc.model.RectangularParallelepiped(
                xmin, xmax, ymin, ymax, zmin, zmax
            )
        elif surface_type == "vessel":
            if len(coefficients) == 7:
                x0, y0, r, zmin, zmax, hbottom, htop = coefficients
            else:
                x0, y0, r, zmin, zmax, hbottom = coefficients
                htop = hbottom
            surface = openmc.model.Vessel(r, zmin, zmax, hbottom, htop, (x0, y0))
        else:
            raise NotImplementedError(
                f"Surface type '{surface_type}' not yet supported in IR conversion."
            )

        if boundary_type and hasattr(surface, "boundary_type"):
            try:
                surface.boundary_type = boundary_type
            except Exception:
                pass

        openmc_surfaces[name] = surface

    return openmc_surfaces


def normalize_region_tokens(tokens: list[str]) -> list[str]:
    normalized = list(tokens)
    for index in range(len(normalized) - 1, 0, -1):
        if normalized[index] == "-":
            normalized[index + 1] = f"-{normalized[index + 1]}"
            del normalized[index]
    return normalized


def cells_from_ir(
    cells_ir: Dict[str, Dict[str, Any]],
    surfaces: Dict[str, openmc.Surface],
    materials: Dict[str, openmc.Material],
    fills: Dict[str, Any] | None = None,
    universes: Dict[str, openmc.Universe] | None = None,
) -> tuple[Dict[str, openmc.Cell], Dict[str, openmc.Universe], set[openmc.Cell]]:
    """Convert cell IR records into OpenMC Cells and Universes."""
    fills = fills or {}
    cells: Dict[str, openmc.Cell] = {}
    universes = universes or {}
    outside_cells: set[openmc.Cell] = set()

    starting_id = openmc.Surface.next_id
    name_to_index: Dict[str, int] = {}
    index_to_surface: Dict[int, openmc.Surface] = {}
    for name, surface in surfaces.items():
        if not name.isnumeric():
            name_to_index[name] = index = starting_id
            starting_id += 1
        else:
            index = int(name)
        index_to_surface[index] = surface

    for name, record in cells_ir.items():
        universe_name = record.get("universe")
        if not universe_name:
            continue

        cell = openmc.Cell(name=name)

        if universe_name not in universes:
            uid = int(universe_name) if universe_name.isnumeric() else None
            universes[universe_name] = openmc.Universe(
                universe_id=uid, name=universe_name
            )
        universes[universe_name].add_cell(cell)

        fill_mode = record.get("fill_mode")
        region_tokens: List[str] = []
        if fill_mode == "fill":
            fill_name = record.get("fill")
            if not fill_name:
                raise ValueError(f"Cell '{name}' missing fill target.")
            fill = fills.get(fill_name)
            if fill is None:
                raise ValueError(
                    f"Cell '{name}' is filled with unknown universe/lattice '{fill_name}'."
                )
            cell.fill = fill
        elif fill_mode == "void":
            pass
        elif fill_mode == "outside":
            outside_cells.add(cell)
        else:
            material_name = record.get("material")
            if material_name:
                if material_name not in materials:
                    raise ValueError(
                        f"Cell '{name}' references missing material '{material_name}'."
                    )
                cell.fill = materials[material_name]

        region_raw = record.get("region", "")
        if region_raw:
            region_tokens = region_raw.split()

        if region_tokens:
            region_tokens = normalize_region_tokens(region_tokens)
            expression = " ".join(region_tokens)
            for surf_name, index in sorted(
                name_to_index.items(), key=lambda item: len(item[0]), reverse=True
            ):
                expression = expression.replace(surf_name, str(index))

            try:
                cell.region = openmc.Region.from_expression(
                    expression=expression, surfaces=index_to_surface
                )
            except Exception as exc:
                raise ValueError(
                    f"Failed to convert cell definition for '{name}': {expression}"
                ) from exc

        cells[name] = cell

    return cells, universes, outside_cells


def get_or_create_universe(
    name: str, universes: Dict[str, openmc.Universe]
) -> openmc.Universe:
    if name in universes:
        return universes[name]
    uid = int(name) if name.isnumeric() else None
    universes[name] = openmc.Universe(universe_id=uid, name=name)
    return universes[name]


def pins_from_ir(
    pins_ir: Dict[str, Dict[str, Any]],
    materials: Dict[str, openmc.Material],
    fills: Dict[str, Any] | None = None,
) -> Dict[str, openmc.Universe]:
    """Convert pin IR records into OpenMC Universes."""
    fills = fills or {}
    universes: Dict[str, openmc.Universe] = {}

    for name, record in pins_ir.items():
        layers = record.get("layers") or []
        if not layers:
            continue

        fills_list: List[Any] = []
        surfaces: List[openmc.Surface] = []
        for layer in layers:
            fill_mode = layer.get("fill_mode", "material")
            fill_name = layer.get("fill")

            if fill_mode == "void":
                fills_list.append(None)
            elif fill_mode == "fill":
                if not fill_name:
                    raise ValueError(f"Pin '{name}' missing layer fill.")
                fill = fills.get(fill_name)
                if fill is None:
                    raise ValueError(
                        f"Pin '{name}' references missing universe '{fill_name}'."
                    )
                fills_list.append(fill)
            else:
                if not fill_name:
                    raise ValueError(f"Pin '{name}' missing layer fill.")
                if fill_name in materials:
                    fills_list.append(materials[fill_name])
                elif fill_name in fills:
                    fills_list.append(fills[fill_name])
                else:
                    raise ValueError(
                        f"Pin '{name}' references missing material/universe '{fill_name}'."
                    )

            radius = layer.get("radius")
            if radius is not None:
                surfaces.append(openmc.ZCylinder(r=float(radius)))

        if surfaces:
            universes[name] = openmc.model.pin(surfaces, fills_list, name=name)
        else:
            cell = openmc.Cell(fill=fills_list[0])
            universes[name] = openmc.Universe(cells=[cell])

    return universes


def build_rectangular_universes(
    names: Iterable[str],
    nx: int,
    ny: int,
    universes: Dict[str, openmc.Universe],
) -> List[List[openmc.Universe]]:
    entries = list(names)
    if len(entries) != nx * ny:
        raise ValueError(
            f"Expected {nx * ny} lattice entries, got {len(entries)}."
        )
    rows = [entries[i * nx : (i + 1) * nx] for i in range(ny)]
    rows = rows[::-1]
    return [[get_or_create_universe(name, universes) for name in row] for row in rows]


def build_rectangular_name_grid(
    names: Iterable[str], nx: int, ny: int, reverse: bool = True
) -> List[List[str]]:
    entries = list(names)
    if len(entries) != nx * ny:
        raise ValueError(
            f"Expected {nx * ny} lattice entries, got {len(entries)}."
        )
    rows = [entries[i * nx : (i + 1) * nx] for i in range(ny)]
    return rows[::-1] if reverse else rows


def choose_outer_name(grid: List[List[str]]) -> str:
    border: List[str] = []
    if not grid:
        raise ValueError("Empty lattice grid.")
    border.extend(grid[0])
    border.extend(grid[-1])
    border.extend(row[0] for row in grid)
    border.extend(row[-1] for row in grid)
    return Counter(border).most_common(1)[0][0]


def rect_to_axial(
    row: int,
    col: int,
    center_row: int,
    center_col: int,
    orientation: str,
) -> Tuple[int, int]:
    """Map Serpent's rectangular lattice indices to axial hex coordinates.

    Serpent hex lattices (type 2/3) list universes on an nx-by-ny rectangle.
    That rectangle represents a rhombus in axial space; we anchor the rectangle
    center at (0, 0) and use hex distance to determine the ring index.
    """
    if orientation == "x":
        q = col - center_col
        r = row - center_row
    else:
        q = row - center_row
        r = col - center_col
    return q, r


def axial_distance(q: int, r: int) -> int:
    s = -q - r
    return max(abs(q), abs(r), abs(s))


def hex_universe_index(
    x: int, a: int, num_rings: int, orientation: str
) -> Tuple[int, int]:
    z = -a - x
    g = max(abs(x), abs(a), abs(z))
    i_ring = num_rings - 1 - g

    if g == 0:
        return i_ring, 0

    if x >= 0:
        if a >= 0:
            i_within = x
        else:
            i_within = 2 * g + z
    else:
        if a <= 0:
            i_within = 3 * g - x
        else:
            i_within = 5 * g - z

    if orientation == "x":
        i_within = (i_within + 5 * g) % (6 * g)

    return i_ring, i_within


def hex_lattice_from_grid(
    name: str,
    grid: List[List[str]],
    pitch: float,
    center: Tuple[float, float],
    orientation: str,
    universes: Dict[str, openmc.Universe],
    lattice_id: int | None,
) -> openmc.HexLattice:
    if not grid or not grid[0]:
        raise ValueError(f"Empty lattice grid for '{name}'.")

    ny = len(grid)
    nx = len(grid[0])
    center_row = ny // 2
    center_col = nx // 2
    # Map each (row, col) in the Serpent lattice listing to axial (q, r)
    # coordinates centered at (0, 0). Rings are based on hex distance.
    centered_map: Dict[Tuple[int, int], str] = {}
    for row_idx, row in enumerate(grid):
        if len(row) != nx:
            raise ValueError(f"Inconsistent lattice row length for '{name}'.")
        for col_idx, entry in enumerate(row):
            q, r = rect_to_axial(
                row_idx, col_idx, center_row, center_col, orientation
            )
            centered_map[(q, r)] = entry

    if not centered_map:
        raise ValueError(f"Empty lattice grid for '{name}'.")

    num_rings = (max(nx, ny) + 1) // 2

    outer_name = choose_outer_name(grid)
    outer_universe = get_or_create_universe(outer_name, universes)

    rings: List[List[openmc.Universe]] = []
    for ring_index in range(num_rings):
        ring_radius = num_rings - 1 - ring_index
        ring_len = max(6 * ring_radius, 1)
        rings.append([outer_universe] * ring_len)

    for (x, a), entry_name in centered_map.items():
        g = axial_distance(x, a)
        if g >= num_rings:
            continue
        i_ring, i_within = hex_universe_index(x, a, num_rings, orientation)
        rings[i_ring][i_within] = get_or_create_universe(entry_name, universes)

    lattice = openmc.HexLattice(lattice_id=lattice_id, name=name)
    lattice.orientation = orientation
    lattice.center = center
    lattice.pitch = [pitch]
    lattice.universes = rings
    lattice.outer = outer_universe
    return lattice


def apply_rectangular_symmetry(
    grid: List[List[str]], usym: Dict[str, Any]
) -> List[List[str]]:
    angle = usym.get("angle")
    if angle is None:
        return grid

    if abs(angle - 45.0) < 1.0e-6:
        transforms = (
            lambda dx, dy: (dx, dy),
            lambda dx, dy: (-dx, dy),
            lambda dx, dy: (dx, -dy),
            lambda dx, dy: (-dx, -dy),
            lambda dx, dy: (dy, dx),
            lambda dx, dy: (-dy, dx),
            lambda dx, dy: (dy, -dx),
            lambda dx, dy: (-dy, -dx),
        )
    elif abs(angle - 90.0) < 1.0e-6:
        transforms = (
            lambda dx, dy: (dx, dy),
            lambda dx, dy: (-dy, dx),
            lambda dx, dy: (-dx, -dy),
            lambda dx, dy: (dy, -dx),
        )
    elif abs(angle - 180.0) < 1.0e-6:
        transforms = (
            lambda dx, dy: (dx, dy),
            lambda dx, dy: (-dx, -dy),
        )
    else:
        raise NotImplementedError(
            f"Unsupported usym angle {angle} for rectangular lattice expansion."
        )

    ny = len(grid)
    nx = len(grid[0]) if ny else 0
    center_row = (ny - 1) / 2.0
    center_col = (nx - 1) / 2.0

    counts = Counter(entry for row in grid for entry in row)
    background = counts.most_common(1)[0][0]

    expanded = [list(row) for row in grid]

    for row_idx, row in enumerate(grid):
        for col_idx, entry in enumerate(row):
            if entry == background:
                continue
            dx = col_idx - center_col
            dy = row_idx - center_row
            for transform in transforms:
                tdx, tdy = transform(dx, dy)
                r = tdy + center_row
                c = tdx + center_col
                r_idx = int(round(r))
                c_idx = int(round(c))
                if abs(r - r_idx) > 1.0e-6 or abs(c - c_idx) > 1.0e-6:
                    continue
                if not (0 <= r_idx < ny and 0 <= c_idx < nx):
                    continue
                current = expanded[r_idx][c_idx]
                if current == background or current == entry:
                    expanded[r_idx][c_idx] = entry
                else:
                    raise ValueError(
                        "Symmetry expansion conflict at "
                        f"({r_idx}, {c_idx}): {current} vs {entry}"
                    )

    return expanded


def lattices_from_ir(
    lattices_ir: Dict[str, Dict[str, Any]],
    universes: Dict[str, openmc.Universe],
    settings: Dict[str, Any] | None = None,
) -> Dict[str, openmc.Lattice]:
    """Convert lattice IR records into OpenMC Lattices."""
    lattices: Dict[str, openmc.Lattice] = {}
    usym_map: Dict[str, Dict[str, Any]] = {}
    if settings:
        for record in settings.get("usym", []):
            universe = record.get("universe")
            if universe:
                usym_map[universe] = record

    for name, record in lattices_ir.items():
        lattice_type = record.get("lat_type_int")
        if lattice_type is None:
            lat_type_raw = record.get("lat_type")
            try:
                lattice_type = int(lat_type_raw) if lat_type_raw is not None else None
            except ValueError:
                lattice_type = None
        if lattice_type is None:
            continue

        lattice_id = int(name) if name.isnumeric() else None

        if lattice_type in (1, 2, 3, 14):
            try:
                x0 = float(record["x0"])
                y0 = float(record["y0"])
                nx = int(record["nx"])
                ny = int(record["ny"])
                pitch = float(record["pitch"])
                entry_names = list(record["entries"])
            except KeyError as exc:
                raise ValueError(
                    f"Lattice '{name}' missing required fields for type {lattice_type}."
                ) from exc
            name_grid = build_rectangular_name_grid(entry_names, nx, ny, reverse=False)
            usym_record = usym_map.get(name)
            if usym_record and lattice_type == 1:
                name_grid = apply_rectangular_symmetry(name_grid, usym_record)
                entry_names = [entry for row in name_grid for entry in row]

            if lattice_type == 1:
                lattice_universes = build_rectangular_universes(
                    entry_names, nx, ny, universes
                )
                lattice = openmc.RectLattice(lattice_id=lattice_id, name=name)
                lattice.lower_left = (x0 - (nx / 2) * pitch, y0 - (ny / 2) * pitch)
                lattice.pitch = (pitch, pitch)
                lattice.universes = lattice_universes
            elif lattice_type in (2, 3):
                orientation = "x" if lattice_type == 2 else "y"
                lattice = hex_lattice_from_grid(
                    name,
                    name_grid,
                    pitch,
                    center=(x0, y0),
                    orientation=orientation,
                    universes=universes,
                    lattice_id=lattice_id,
                )
            else:
                raise NotImplementedError(
                    "X-type triangular lattice not yet supported."
                )

        elif lattice_type in (6, 7, 8):
            try:
                x0 = float(record["x0"])
                y0 = float(record["y0"])
                pitch = float(record["pitch"])
                universe_name = record["entries"][0]
            except (KeyError, IndexError) as exc:
                raise ValueError(
                    f"Lattice '{name}' missing required fields for type {lattice_type}."
                ) from exc
            universe = get_or_create_universe(universe_name, universes)

            if lattice_type == 6:
                lattice = openmc.RectLattice(lattice_id=lattice_id, name=name)
                lattice.lower_left = (-(x0 + pitch / 2), -(y0 + pitch / 2))
                lattice.pitch = (pitch, pitch)
                lattice.universes = [[universe]]
                lattice.outer = universe
            elif lattice_type in (7, 8):
                lattice = openmc.HexLattice(lattice_id=lattice_id, name=name)
                lattice.orientation = "x" if lattice_type == 7 else "y"
                lattice.center = (x0, y0)
                lattice.pitch = [pitch]
                lattice.universes = [[universe]]
                lattice.outer = universe

        elif lattice_type == 9:
            try:
                x0 = float(record["x0"])
                y0 = float(record["y0"])
                n = int(record["n"])
                z_values = [float(value) for value in record["z_values"]]
                uni_names = list(record["stack_universes"])
            except KeyError as exc:
                raise ValueError(
                    f"Lattice '{name}' missing required fields for type 9."
                ) from exc
            if len(z_values) != n:
                raise ValueError(f"Expected {n} z-values for lattice '{name}'.")
            stack_universes = [
                get_or_create_universe(unit, universes) for unit in uni_names
            ]
            lattice = vertical_stack(z_values, stack_universes, x0, y0)
        else:
            raise NotImplementedError(
                f"Lattice type '{lattice_type}' not yet supported."
            )

        lattices[name] = lattice
        universes[name] = lattice

    return lattices


def geometry_components_from_ir(
    geom_ir: Dict[str, Any], materials: Dict[str, openmc.Material]
) -> Dict[str, Any]:
    """Convert geometry IR records into OpenMC objects and lookup maps."""
    surfaces = surfaces_from_ir(geom_ir["surfaces"])
    pins = pins_from_ir(geom_ir["pins"], materials)

    universes: Dict[str, openmc.Universe] = dict(pins)
    lattices = lattices_from_ir(
        geom_ir["lattices"], universes, settings=geom_ir.get("settings")
    )

    fills: Dict[str, Any] = dict(universes)
    fills.update(lattices)

    cells, cell_universes, outside_cells = cells_from_ir(
        geom_ir["cells"], surfaces, materials, fills, universes=universes
    )

    return {
        "surfaces": surfaces,
        "pins": pins,
        "lattices": lattices,
        "cells": cells,
        "universes": cell_universes,
        "outside_cells": outside_cells,
    }
