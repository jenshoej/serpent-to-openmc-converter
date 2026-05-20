# SPX `run-spx.py` Rundown

This document explains what `run-spx.py` creates on top of the Serpent-to-OpenMC conversion, and what parts are good candidates for moving into the converter.

## What `run-spx.py` creates

## 1. Base converted model (already converter-level)
- Builds OpenMC materials + geometry from Serpent `mat`, `surf`, `cell`, `pin`, and `lat` cards via `build_openmc_components(...)`.
- Creates `openmc.Model(geometry=..., materials=...)`.

## 2. Boundary handling
- Finds outside cells and sets their outer surfaces to `vacuum` when they are unset/transmission.
- This is currently a post-conversion override in the script.

## 3. Geometry extents and plotting bounds
- Computes finite bounding boxes for source and plotting.
- Falls back through multiple strategies when geometry bounding boxes are non-finite.
- Infers an XY plot z-plane from Serpent type-9 stack lattices (SPX naming convention, `u_02...` as core section).

## 4. Plots
- Creates two OpenMC plots:
  - `spx_xy` (basis `xy`)
  - `spx_xz` (basis `xz`)
- Uses material coloring and writes images with `model.plot_geometry(...)`.

## 5. Tallies
- Adds:
  - `global_power` (`kappa-fission`, `fission`)
  - `global_reaction_rates` (`kappa-fission`, `fission`, `nu-fission`, `absorption`, `flux`)
  - `kinetics_ifp` (`ifp-time-numerator`, `ifp-beta-numerator`, `ifp-denominator`)
  - `xy_power_map` on a 20x20x1 mesh (`kappa-fission`)
  - `xy_reaction_map` on same mesh (reaction scores above)
  - `xz_power_map` on a 20x1x20 mesh (`kappa-fission`)
  - `xz_reaction_map` on same mesh (reaction scores above)
  - `cell_power_<id>` tallies for up to 6 sampled fissionable cells
  - `subassembly_power_cells` tally plus `subassembly_power_map.json` metadata for
    reconstructing a full assembly-wise power map from the `l_core` lattice

Notes:
- `k_eff` is not a tally; it is produced by eigenvalue mode by default.
- Prompt lifetime and `beta_eff` are obtained from the IFP tally components in post-processing.

## 6. Settings
- Sets OpenMC eigenvalue settings:
  - `run_mode = eigenvalue`
  - particles / batches / inactive
  - `ifp_n_generation`
  - temperature method/range
  - box source with `only_fissionable=True`

## 7. Outputs and run
- Exports XML (`materials.xml`, `geometry.xml`, `settings.xml`, `tallies.xml`, `plots.xml`, `model.xml`).
- Writes `subassembly_power_map.json` when the assembly-wise power map metadata is available.
- Runs OpenMC.
- Finds latest statepoint and prints `k-effective`.

## OpenMC data library setup

`run-spx.py` requires `OPENMC_CROSS_SECTIONS` to point at a readable OpenMC
`cross_sections.xml`.

In this workspace, you can build a JEFF-3.1.1 HDF5 library from the processed
ACE files in `/home/christiansen/JEFF311-Processed-ACE/N-ACE` with:

```bash
micromamba run -n openmc python /data/workspace/ACT/ACT-26-02-SERPENT-TO-OPENMC/openmc_serpent_adapter/SPX_case/convert_jeff311_ace.py
```

That creates the library at:

```text
/data/workspace/ACT/ACT-26-02-SERPENT-TO-OPENMC/openmc_serpent_adapter/SPX_case/openmc_data/jeff311-hdf5
```

To activate it for the current shell:

```bash
source /data/workspace/ACT/ACT-26-02-SERPENT-TO-OPENMC/openmc_serpent_adapter/SPX_case/use_jeff311.sh
```

That helper now sets both:

- `OPENMC_CROSS_SECTIONS` to the JEFF-3.1.1 HDF5 library
- `SPX_OUTPUT_DIR` to `SPX_case/jeff3.1.1` so the new run stays separate from archived results

Then run the case:

```bash
micromamba run -n openmc python /data/workspace/ACT/ACT-26-02-SERPENT-TO-OPENMC/openmc_serpent_adapter/SPX_case/run-spx.py
```

`conversion_report.json` in the output directory lists the JEFF311 duplicate
table-name collisions that OpenMC cannot distinguish directly from the ACE IDs.

---

## Cross-check against `SPX` input

Top-level cards in `SPX` are:
- `mat`, `surf`, `cell`, `pin`, `lat`, `set`, `det`, `plot`, `ene`

Current converter path covers:
- `mat`, `surf`, `cell`, `pin`, `lat` (core geometry/material conversion)
- `set usym` only (when present)

Not converted from this `SPX` file:
- `det` cards (67 detector definitions)
- `plot` cards (22 plot definitions)
- Most `set` cards (`pop`, `power`, `seed`, `bc`, `egrid`, `micro`, `nfg`, etc.)
- `ene` energy grid card

`run-spx.py` currently replaces these with OpenMC-native setup (custom tallies/settings/plots).

---

## What can be added to the converter fairly easily

## High value, low effort
1. `set pop` -> OpenMC settings mapping
- Serpent: `set pop NPG NGEN NSKIP`
- OpenMC:
  - `particles = NPG`
  - `inactive = NSKIP`
  - `batches = NGEN + NSKIP`

2. `set seed` -> OpenMC RNG seed
- Map directly to `settings.seed`.

3. `set bc` -> geometry boundary policy
- Add a converter option to apply `vacuum`/periodic behavior from Serpent BC intent, instead of hardcoding in case scripts.

## Medium effort (worth doing next)
4. `det` to OpenMC tally translation (subset first)
- Start with detector patterns used in SPX (`du`, `dr`, `dz`, lattice/universe-based power/fission profiles).
- This removes most case-specific tally code.

5. `plot` card translation (subset first)
- Translate common SPX `plot` forms into OpenMC `Plot` objects.

## Lower priority / limited direct mapping
6. `set power`
- Not a direct OpenMC eigenvalue setting field.
- Useful for normalization/depletion workflows, but needs a policy (post-processing normalization vs depletion normalization).

7. `set micro` + `ene`
- These are Serpent multigroup/generation inputs; no direct one-to-one in this CE OpenMC conversion path.

---

## Practical split of responsibilities

Recommended converter responsibility:
- Convert geometry/materials.
- Convert straightforward global settings (`pop`, `seed`, selected `bc` behavior).
- Optionally convert a supported detector/plot subset.

Recommended case-script responsibility:
- SPX-specific heuristics (e.g., type-9 section inference by universe naming).
- Analysis-specific tallies and plotting choices not explicitly defined by converted Serpent cards.
