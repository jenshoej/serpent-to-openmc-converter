
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

import openmc

from .ir_openmc import geometry_components_from_ir, materials_from_ir
from .serpent_geom_ir import parse_geometry_ir
from .serpent_mat_ir import parse_materials_to_ir, parse_therm_cards
from .serpent_preprocess import load_serpent_lines


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
    materials, geom_components, root = build_openmc_components(
        input_path, root_universe=root_universe
    )

    model = openmc.Model(geometry=openmc.Geometry(root))
    model.materials = openmc.Materials(materials.values())
    return model


def build_openmc_components(
    input_path: Path, root_universe: str | None = None
) -> tuple[Dict[str, openmc.Material], Dict[str, object], openmc.Universe]:
    """Return OpenMC materials/geometry objects for manual editing."""
    lines = load_serpent_lines(input_path)

    therm = parse_therm_cards(lines)
    materials_ir = parse_materials_to_ir(lines, therm, preprocess=False)
    materials = materials_from_ir(materials_ir)

    geom_ir = parse_geometry_ir(lines)
    geom_components = geometry_components_from_ir(geom_ir, materials)
    root = select_root_universe(
        root_universe, geom_components["universes"], geom_components["outside_cells"]
    )

    return materials, geom_components, root


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
