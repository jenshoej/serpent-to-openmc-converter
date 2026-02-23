# Serpent Conversion Tools for OpenMC

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

This repository provides tools for parsing/converting Serpent models to OpenMC
classes and/or XML files.

## What This Produces

Given a Serpent input file, the converter builds an OpenMC geometry/materials
model from supported cards and can produce:

- An in-memory `openmc.Model` object (Python API).
- A `model.xml` file, written with `openmc.Model.export_to_model_xml(...)`.
- Optional geometry plot images if you call `plot_model(...)`.

The converter currently targets material and geometry translation. It does not
automatically create OpenMC source/tally/settings definitions equivalent to
Serpent `src`/`det` input.

## How To Use It

### 1. Install in editable mode

From the repository root:

```bash
pip install -e .
```

### 2. Convert and extend in Python

```python
from pathlib import Path
import openmc
from src.serpent_to_openmc import build_openmc_components

materials, geom_components, root = build_openmc_components(Path("path/to/input"))
model = openmc.Model(geometry=openmc.Geometry(root))
model.materials = openmc.Materials(materials.values())

# Add OpenMC-native features (example: settings/tallies)
model.settings = openmc.Settings()
model.settings.batches = 50
model.settings.inactive = 10
model.settings.particles = 10000

tally = openmc.Tally(name="flux")
tally.filters = [openmc.CellFilter(list(geom_components["cells"].values()))]
tally.scores = ["flux"]
model.tallies = openmc.Tallies([tally])

model.export_to_model_xml(path="model.xml")
```

`build_openmc_components(...)` returns:

- `materials`: `dict[str, openmc.Material]`
- `geom_components`: dict containing `surfaces`, `pins`, `lattices`, `cells`,
  `universes`, and `outside_cells`
- `root`: selected root `openmc.Universe`

`build_openmc_model(...)` is also available if you just want a ready-to-export
model before adding your own settings/tallies.

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

- Only `set usym` is parsed. `set root` and `set bc` are ignored (manual selection).
