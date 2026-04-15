# P3 Offline Post-Processor

Computes diagnostic ice microphysics quantities from P3 prognostic variables using lookup table (LUT) interpolation, without running the full model.

## Motivation

When running P3 in a full-scale NWP model (e.g., GEM), only the **prognostic variables** are output:
- `Qitot` — total ice mass mixing ratio
- `Qirim` — rime ice mass mixing ratio  
- `Nitot` — total ice number concentration
- `Birim` — rime volume (bulk density = Qirim / Birim)
- `Qiliq` — liquid water in ice (optional, if liqFrac enabled)
- `Zitot` — reflectivity moment (optional, if 3-moment ice enabled)

The **diagnostic quantities** (e.g., bulk ice density, mean-mass diameter, effective radius, and PSD parameters) are normally computed internally by `p3_main`. This post-processor replicates that computation offline, so users with only prognostic output can recover the diagnostics without re-running the model.

## Quick Start

From `./kin1d/`:

```bash
python3 postproc/p3_postproc.py
```

This reads:
- `out_p3.dat` (model output from cld1d)
- `lookup_tables/p3_lookupTable_*.dat-v*` (unzipped LUT files)
- `p3_config.txt` (P3 configuration flags)

And writes:
- `out_pp.dat` — same format as `out_p3.dat`, with six extra columns appended per ice category:
  - `pp_rhoi` — post-processor bulk ice density [kg m-3]
  - `pp_di` — post-processor mean-mass diameter [m]
  - `pp_effi` — post-processor effective radius [m]
  - `pp_vmi` — post-processor mass-weighted fall speed [m s-1], using `rhofaci = 1`
  - `pp_lami` — post-processor PSD slope parameter lambda_i [m-1]
  - `pp_mui` — post-processor PSD shape parameter mu_i [-]

The original model `diag_rhoi` and `diag_di` columns are preserved, allowing direct comparison. If the `cld1d` output also includes the optional verification fields for `diag_effi`, `diag_vmi`, `diag_lami`, and `diag_mui`, the script will compare those as well.

## Command-Line Options

```bash
python3 postproc/p3_postproc.py \
  --input FILE           # Input data file (default: out_p3.dat)
  --lut-dir DIR          # Directory with unzipped LUT files (default: lookup_tables)
  --config FILE          # P3 config file (default: p3_config.txt)
  --output FILE          # Output file (default: out_pp.dat)
```

## Supported Configurations

**Currently implemented:**
- `nCat = 1..3`
- `trplMomIce = .false.` and `.true.`
- `liqFrac = .false.` and `.true.`
- Input from `out_p3.dat` via CLI
- Direct NumPy-array input via the `compute()` Python API

**Not yet implemented:**
- netCDF input
- GEM FST/RPN input
- Additional diagnostics that require extra state beyond the current LUT lookups, such as `diag_ZET`

## Verification / Development

The script includes built-in statistics to compare post-processor output against the model's original diagnostics. After running, it prints:
- Number of active ice levels processed per category
- Max and mean absolute differences
- Max and median relative errors

**Expected accuracy:** Median relative error is typically at the ASCII precision floor, about `0.01-0.03%` for both `diag_rhoi` and `diag_di`, across all supported configurations. Remaining max-error outliers occur at a small number of 3-moment levels where `apply_mui_bounds_to_zi` in `p3_main` modified `Zitot` before file output.

## Implementation Details

The post-processor replicates the key Fortran routines from `microphy_p3.f90`:

1. **`calc_bulkRhoRime`** — Compute bulk rime density from prognostic variables:
   ```
   rho_rime = Qirim / Birim  (with bounds checking)
   ```

2. **`find_lookupTable_indices_1a`** — Map prognostic variables to LUT indices via interpolation formulas:
   ```
   dumi  = f(log10(Qitot/Nitot))     # particle size index
   dumii = f(Qirim/(Qitot-Qiliq))    # rime fraction index
   dumll = f(Qiliq/Qitot)            # liquid fraction index
   dumjj = f(rho_rime)               # density index (non-linear)
   ```

3. **`proc_from_LUT_main2mom`** — 4D (or 3D) linear interpolation in the 2-moment LUT:
   ```
   itab[densize, rimsize, liqsize, isize, tabsize]
   result = interpolate(size, rime, density, liquid_fraction)
   ```

4. **`find_lookupTable_indices_3a`** / **`find_lookupTable_indices_1c`** — additional 3-moment index calculations for `Zi/Qi` and `mu_i`

5. **`proc_from_LUT_3`** and **`proc_from_LUT_main3mom`** — 5D/6D interpolation in the 3-moment LUTs:
   ```
   itab_3mom_mui[zqsize, densize, rimsize, liqsize, isize, 2]
   itab_3mom[zzsize, densize, rimsize, liqsize, isize, tabsize_3mom]
   ```

The current script computes:
- `diag_rhoi`
- `diag_di`
- `diag_effi`
- `diag_vmi` with `rhofaci = 1`
- `diag_lami`
- `diag_mui`

All arithmetic and indexing matches the Fortran code exactly, ensuring reproducibility.

## LUT Files

The script expects decompressed LUT files in the `lookup_tables/` directory. Download from the P3 repository and unzip:

```bash
gunzip lookup_tables/p3_lookupTable_1.dat-v6.9-2momI.gz
gunzip lookup_tables/p3_lookupTable_1.dat-v6.9-3momI.gz
gunzip lookup_tables/p3_lookupTable_2.dat-v6.2.gz
gunzip lookup_tables/p3_lookupTable_3.dat-v1.4.gz
```

Only the tables matching your configuration are required at runtime:
- 2-moment runs need `p3_lookupTable_1.dat-v6.9-2momI`
- 3-moment runs need `p3_lookupTable_1.dat-v6.9-3momI` and `p3_lookupTable_3.dat-v1.4`

## Dependencies

- Python 3.6+
- NumPy

Install NumPy if needed:
```bash
pip3 install numpy
```

## References

- P3 microphysics scheme: Morrison & Milbrandt (2015) and updates
- Lookup table generation and structure: See `./src/create_p3_lookupTable_*.f90`
- Model integration: See `./src/p3/microphy_p3.f90` (p3_init, p3_main)
