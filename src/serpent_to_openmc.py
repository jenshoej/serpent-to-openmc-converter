
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

import openmc

from .ir_openmc import geometry_components_from_ir, materials_from_ir
from .serpent_geom_ir import parse_geometry_ir
from .serpent_mat_ir import parse_materials_to_ir, parse_therm_cards
from .serpent_preprocess import first_word, load_serpent_lines


@dataclass(slots=True)
class SerpentRunSettings:
    """High-level summary of run-specific Serpent `set` cards."""

    particles: int | None = None
    active_generations: int | None = None
    inactive_generations: int | None = None
    suggested_openmc_batches: int | None = None
    seed: int | None = None
    power: float | None = None
    boundary_conditions: tuple[str, ...] = ()
    acelib: str | None = None
    declib: str | None = None
    nfylib: str | None = None
    extras: dict[str, tuple[str, ...]] = field(default_factory=dict)
    raw_set_cards: dict[str, tuple[tuple[str, ...], ...]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "particles": self.particles,
            "active_generations": self.active_generations,
            "inactive_generations": self.inactive_generations,
            "suggested_openmc_batches": self.suggested_openmc_batches,
            "seed": self.seed,
            "power": self.power,
            "boundary_conditions": list(self.boundary_conditions),
            "acelib": self.acelib,
            "declib": self.declib,
            "nfylib": self.nfylib,
            "extras": {key: list(value) for key, value in self.extras.items()},
            "raw_set_cards": {
                key: [list(values) for values in records]
                for key, records in self.raw_set_cards.items()
            },
        }


@dataclass(slots=True)
class ConversionReport:
    """Simple, user-facing summary of what the converter produced."""

    input_path: Path
    root_universe_name: str | None
    material_count: int
    surface_count: int
    cell_count: int
    universe_count: int
    lattice_count: int
    pin_count: int
    outside_cell_names: tuple[str, ...]
    applied_boundary_type: str | None
    boundary_surfaces_changed: int
    run_settings: SerpentRunSettings
    lattice_maps: dict[str, dict[str, Any]] = field(default_factory=dict)
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "input_path": str(self.input_path),
            "root_universe_name": self.root_universe_name,
            "material_count": self.material_count,
            "surface_count": self.surface_count,
            "cell_count": self.cell_count,
            "universe_count": self.universe_count,
            "lattice_count": self.lattice_count,
            "pin_count": self.pin_count,
            "outside_cell_names": list(self.outside_cell_names),
            "applied_boundary_type": self.applied_boundary_type,
            "boundary_surfaces_changed": self.boundary_surfaces_changed,
            "run_settings": self.run_settings.to_dict(),
            "lattice_maps": self.lattice_maps,
            "notes": list(self.notes),
        }


def collect_set_cards(lines: list[str]) -> dict[str, list[list[str]]]:
    """Collect all Serpent `set` cards as token lists keyed by setting name."""
    settings: dict[str, list[list[str]]] = {}
    for line in lines:
        words = line.split()
        if not words or first_word(words) != "set" or len(words) < 2:
            continue
        settings.setdefault(words[1].lower(), []).append(words[2:])
    return settings


def summarize_run_settings_from_lines(lines: list[str]) -> SerpentRunSettings:
    """Summarize run-relevant Serpent `set` cards for user inspection."""
    raw_cards = collect_set_cards(lines)
    known_extra_keys = {
        "opti",
        "egrid",
        "ures",
        "sym",
        "nfg",
        "micro",
    }

    summary = SerpentRunSettings(
        raw_set_cards={
            key: tuple(tuple(values) for values in records)
            for key, records in raw_cards.items()
        }
    )

    pop_records = raw_cards.get("pop", [])
    if pop_records:
        values = pop_records[-1]
        if len(values) >= 3:
            summary.particles = int(values[0])
            summary.active_generations = int(values[1])
            summary.inactive_generations = int(values[2])
            summary.suggested_openmc_batches = (
                summary.active_generations + summary.inactive_generations
            )

    seed_records = raw_cards.get("seed", [])
    if seed_records and seed_records[-1]:
        summary.seed = int(seed_records[-1][0])

    power_records = raw_cards.get("power", [])
    if power_records and power_records[-1]:
        summary.power = float(power_records[-1][0])

    bc_records = raw_cards.get("bc", [])
    if bc_records:
        summary.boundary_conditions = tuple(bc_records[-1])

    for key in ("acelib", "declib", "nfylib"):
        records = raw_cards.get(key, [])
        if records and records[-1]:
            setattr(summary, key, records[-1][0])

    for key in known_extra_keys:
        records = raw_cards.get(key, [])
        if records:
            summary.extras[key] = tuple(records[-1])

    return summary


def summarize_run_settings(input_path: Path) -> SerpentRunSettings:
    """Load a Serpent input file and summarize run-relevant `set` cards."""
    lines = load_serpent_lines(input_path)
    return summarize_run_settings_from_lines(lines)


def _build_openmc_components_from_lines(
    lines: list[str], root_universe: str | None = None
) -> tuple[Dict[str, openmc.Material], Dict[str, object], openmc.Universe]:
    therm = parse_therm_cards(lines)
    materials_ir = parse_materials_to_ir(lines, therm, preprocess=False)
    materials = materials_from_ir(materials_ir)

    geom_ir = parse_geometry_ir(lines)
    geom_components = geometry_components_from_ir(geom_ir, materials)
    root = select_root_universe(
        root_universe, geom_components["universes"], geom_components["outside_cells"]
    )

    return materials, geom_components, root


def build_conversion_report(
    input_path: Path,
    materials: Dict[str, openmc.Material],
    geom_components: Dict[str, object],
    root: openmc.Universe,
    run_settings: SerpentRunSettings,
) -> ConversionReport:
    outside_cells = sorted(
        (cell.name for cell in geom_components["outside_cells"]),
        key=str,
    )
    notes: list[str] = []
    boundary_type = geom_components.get("boundary_type")
    if boundary_type == "vacuum":
        notes.append("Applied Serpent `set bc 1` as OpenMC vacuum boundaries.")
    elif run_settings.boundary_conditions:
        notes.append(
            "Serpent boundary conditions were detected, but only `set bc 1` "
            "is currently applied automatically."
        )
    if any(
        value is not None
        for value in (
            run_settings.particles,
            run_settings.active_generations,
            run_settings.inactive_generations,
            run_settings.seed,
            run_settings.power,
        )
    ):
        notes.append(
            "Run-relevant Serpent `set` cards were summarized in "
            "`report.run_settings` for manual OpenMC settings setup."
        )
    if geom_components.get("lattice_maps"):
        notes.append(
            "Serpent lattice-entry maps are available in `report.lattice_maps` "
            "and `geom_components['lattice_maps']` for assembly-wise postprocessing."
        )

    return ConversionReport(
        input_path=input_path,
        root_universe_name=root.name,
        material_count=len(materials),
        surface_count=len(geom_components["surfaces"]),
        cell_count=len(geom_components["cells"]),
        universe_count=len(geom_components["universes"]),
        lattice_count=len(geom_components["lattices"]),
        pin_count=len(geom_components["pins"]),
        outside_cell_names=tuple(outside_cells),
        applied_boundary_type=boundary_type,
        boundary_surfaces_changed=int(
            geom_components.get("boundary_surfaces_changed", 0)
        ),
        run_settings=run_settings,
        lattice_maps=geom_components.get("lattice_maps", {}),
        notes=tuple(notes),
    )


def select_root_universe(
    preferred: str | None,
    universes: Dict[str, openmc.Universe],
    outside_cells: Iterable[openmc.Cell],
) -> openmc.Universe:
    if preferred:
        if preferred not in universes:
            raise ValueError(f"Root universe '{preferred}' not found.")
        return universes[preferred]

    if "0" in universes:
        return universes["0"]

    if len(universes) == 1:
        return next(iter(universes.values()))

    outside_set = set(outside_cells)
    for universe in universes.values():
        if any(cell in outside_set for cell in universe.cells.values()):
            return universe

    raise ValueError(
        "Unable to determine root universe. Pass root_universe to specify one."
    )


def build_openmc_model(
    input_path: Path, root_universe: str | None = None
) -> openmc.Model:
    model, _ = build_model(input_path, root_universe=root_universe)
    return model


def build_model(
    input_path: Path, root_universe: str | None = None
) -> tuple[openmc.Model, ConversionReport]:
    """Build a base OpenMC model plus a simple conversion report.

    This is the most convenient entry point when you want a ready-to-edit
    OpenMC model and an at-a-glance summary of what was converted from Serpent.
    """
    lines = load_serpent_lines(input_path)
    materials, geom_components, root = _build_openmc_components_from_lines(
        lines, root_universe=root_universe
    )
    run_settings = summarize_run_settings_from_lines(lines)

    model = openmc.Model(geometry=openmc.Geometry(root))
    model.materials = openmc.Materials(materials.values())
    report = build_conversion_report(
        input_path,
        materials,
        geom_components,
        root,
        run_settings,
    )
    return model, report


def build_openmc_components(
    input_path: Path, root_universe: str | None = None
) -> tuple[Dict[str, openmc.Material], Dict[str, object], openmc.Universe]:
    """Return OpenMC materials/geometry objects for manual editing."""
    lines = load_serpent_lines(input_path)
    return _build_openmc_components_from_lines(lines, root_universe=root_universe)


def plot_model(
    model: openmc.Model,
    basis: str,
    origin: Tuple[float, float, float],
    width: Tuple[float, float],
    pixels: Tuple[int, int],
    color_by: str,
    filename: str,
    output_dir: Path,
    openmc_exec: str,
) -> None:
    plot = openmc.Plot()
    plot.basis = basis
    plot.origin = origin
    plot.width = width
    plot.pixels = pixels
    plot.color_by = color_by
    plot.filename = filename
    model.plots = openmc.Plots([plot])

    model.plot_geometry(
        cwd=output_dir,
        openmc_exec=openmc_exec,
        export_model_xml=False,
        directory=output_dir,
    )
