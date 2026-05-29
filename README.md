# Serpent to OpenMC

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

This repository provides a Python tool for parsing/converting Serpent input files to OpenMC. It is a continuation of an unfinished version from the openmc-dev team  (https://github.com/openmc-dev/openmc_serpent_adapter).

Given a Serpent input file, the converter builds an OpenMC
model containing the geometry, materials and a handful of other features from supported cards (not all cards are supported; see below).

In other words, it is not a complete converter translating every Serpent input card, like more run-specific `src` and `det` cards.

The desired OpenMC run settings should be created manually (through the Python API) after the converting the Serpent input, before exporting and running the final OpenMC model.


## How To Use It

The code includes four user functions:

`build_openmc_model()`
- returns: model (OpenMC Class), report (`ConversionReport`)
- The model contains the materials, geometry and supported features of the OpenMC Model.
- The report serves as a description of what the converted OpenMC model contains.
- This function is the default main converter.

`build_openmc_components()`
- returns: materials, geometry, stuff
- Used if the

`summarize_run_settings()`
- returns: a report containing the run settings from the Serpent file. Can be useful to determine the tallies and run settings for the OpenMC model.

`plot_model()`
- Optional function to generate plot(s) of the OpenMC model

## Demonstration

This code snippet demonstrates using the main converter function `build_openmc_model()`. The returned `report` gives a compact description of what was converted.

```python
from pathlib import Path
import openmc
from src.serpent_to_openmc import build_openmc_model

# Run build_openmc_model
model, report = build_openmc_model(Path("path/to/Serpent_input"))

# Inspect what the conversion produced
print(report)
print(report.describe())
print(report.to_dict())
print(report.run_settings)

# Add desired OpenMC features (settings, tallies etc.)
model.settings = openmc.Settings()
model.settings.batches = 50
model.settings.inactive = 10
model.settings.particles = 10000

tally = openmc.Tally(name="flux")
tally.filters = [openmc.CellFilter(list(model.geometry.get_all_cells().values()))]
tally.scores = ["flux"]
model.tallies = openmc.Tallies([tally])

model.export_to_model_xml(path="model.xml")
```

If you prefer direct access to the converted geometry, `build_openmc_components()`
can be used to return materials, geometry-component lookup maps, and the
selected root universe.

Serpent `set` cards are not converted automatically (except a few ones like vacuum boundaries). `summarize_run_settings()` can be used to inspect the Serpent settings if needed in the OpenMC model.

```python
from pathlib import Path
from src.serpent_to_openmc import summarize_run_settings

settings = summarize_run_settings(Path("path/to/input"))
print(settings.particles, settings.inactive_generations, settings.seed)

# A core map can be viewed from the report with
core_map = report.lattice_maps.get("l_core")
if core_map is not None:
    print(core_map["dimension"])
    print(core_map["serpent_grid"][0][0])
```

`build_openmc_components(...)` returns:

- `materials`: `dict[str, openmc.Material]`
- `geom_components`: dict containing `surfaces`, `pins`, `lattices`, `cells`,
  `universes`, `outside_cells`, and `lattice_maps`
- `root`: selected root `openmc.Universe`

`geom_components["lattice_maps"]` preserves Serpent lattice-entry metadata in a
postprocessing-friendly form. For 2D lattices this includes dimensions, pitch,
the Serpent-order `serpent_grid`, and per-position records. Hex lattices also
include axial `(q, r)` coordinates and ring indices.

`build_openmc_model(...)` returns:

- `model`: base `openmc.Model` with converted geometry/materials
- `report`: `ConversionReport` summarizing counts, root universe, outside cells,
  applied boundary condition, parsed run settings, and preserved Serpent
  lattice maps

To describe the converted model after calling `build_openmc_model(...)`, use the returned `report`:

- `print(report)` or `print(report.describe())` for a human-readable summary
- `report.to_dict()` for a structured dictionary view
- `report.run_settings` for the parsed Serpent `set` cards
- `report.lattice_maps` for preserved Serpent lattice-entry metadata

## Logic

The logic of the converter is structured as follows:

Serpent input -> IR layer as Python dicts -> OpenMC equivalent classes

## Known Limitations

The converter currently only handles geometry and material information; source
definition (`src`) and detectors (`det`) should be manually included in the OpenMC model according to it's native syntax (Tallies, Settings)

Supported Serpent cards (IR pipeline):

- `include`, `mat`, `mix`, `therm`
- `surf`, `cell`, `pin`, `lat`
- `set usym` (rectangular symmetry expansion only)

Unsupported/ignored Serpent cards (explicitly no effect):

- `branch`, `casematrix`, `coef`, `datamesh`, `dep`, `det`, `div`, `dtrans`,
  `ene`, `ftrans`, `fun`, `hisv`, `ifc`, `mesh`, `mflow`, `nest`, `particle`,
  `pbed`, `phb`, `plot`, `rep`, `sample`, `sens`, `solid`, `src`, `strans`,
  `thermstoch`, `tme`, `umsh`, `utrans`, `wwgen`, `wwin`

Transformations:

- `trans`, `transa`, `transb`, `transv` are parsed but currently not applied.
- `dtrans`, `ftrans`, `strans`, `utrans` are ignored.

Lattices:

- Supported lattice types: `1`, `2`, `3`, `6`, `7`, `8`, `9`.
- Unsupported lattice types: any other type (including `4`, `11`, `12`, `13`, `14`).

Set cards:

- `set usym` is parsed for rectangular symmetry expansion.
- `set bc 1` is translated to OpenMC `vacuum` on surfaces referenced by Serpent
  `outside` cells.
- `set root` is ignored, and `set bc` modes other than `1` are not yet applied.
