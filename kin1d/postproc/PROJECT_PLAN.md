# P3 Post-Processor — Project Plan

## Goal

Provide an offline tool that computes P3 ice diagnostic quantities from the prognostic
variables output by a model run (GEM, kin1d, or raw arrays), supporting all P3
configurations, without requiring a live model call.

**Prognostic inputs:** `Qitot`, `Qirim`, `Nitot`, `Birim` (always) + `Qiliq` (liqFrac=T) + `Zitot` (3-moment)  
**Diagnostic outputs:** `diag_rhoi`, `diag_di` (initially), expanding to others in Phase 4

---

## Phases

### Phase 1 — Configuration generalization (2-moment ice)

- [x] **1.1** `nCat=1`, `trplMomIce=F`, `liqFrac=F` ← **DONE**
  - Verified: median relative error ~0.03% vs model output
  - Script: `kin1d/postproc/p3_postproc.py`

- [x] **1.2** `nCat=1`, `trplMomIce=F`, `liqFrac=T`
  - Fixed config parser (value ! key format); 4D interpolation path verified

- [x] **1.3** `nCat=2`, `trplMomIce=F`, any liqFrac
  - `pp_rhoi`/`pp_di` now 2D `[n_levels, n_cat]`; output appends `2×n_cat` columns

- [x] **1.4** `nCat=3`, `trplMomIce=F`, any liqFrac
  - Verified; all 3 categories independently produce median rel err ~0.03%

---

### Phase 2 — 3-moment ice  ← **DONE**

- [x] **2.1** Load `itab_3mom` from `p3_lookupTable_1.dat-v6.9-3momI`
  - Shape: `[11, 5, 4, 4, 50, 29]`; loop: `zz→jj→ii→ll→i`; 36 cols/row (5 dummy + 29 + 2)
  - After each `(zz,jj,ii,ll)` block: skip `50×30` collection rows

- [x] **2.2** Load `itab_3mom_mui` from `p3_lookupTable_3.dat-v1.4`
  - Shape: `[80, 5, 4, 4, 50, 2]`; loop: `zq→jj→ii→ll→i`; 9 cols/row; no collection rows

- [x] **2.3** Implement `_find_indices_3a(zi_tot, qi_tot)` → `(dumzq, dum8)`
  - `dum8 = (log10(zi_tot/qi_tot) + 23.) * 3.10347652`

- [x] **2.4** Implement `_find_indices_1c(mu_i)` → `(dumzz, dum6)`
  - `dum6 = mu_i * 0.5 + 1.`

- [x] **2.5** Implement `_interp_3mom_mui()` — 5D interpolation in `itab_3mom_mui`
  - Mirrors `proc_from_LUT_3()` in `microphy_p3.f90`
  - Index 1 → `mu_i`; index 2 → `diag_rhoi`

- [x] **2.6** Implement `_interp_3mom()` — 6D (or 5D if no liqFrac) interpolation in `itab_3mom`
  - Mirrors `proc_from_LUT_main3mom()` in `microphy_p3.f90`
  - Index 11 → `diag_di`

- [x] **2.7** Implement `get_diag_3mom()` top-level call:
  ```
  calc_bulk_rho_rime → _find_indices → _find_indices_3a →
  _interp_3mom_mui(1) [mu_i] → _interp_3mom_mui(2) [rhoi] →
  _find_indices_1c → _interp_3mom(11) [di]
  ```

- [x] **2.8** Route `get_diag()` by `cfg.triple_moment`; pass `Zitot` from input

- [x] **2.9** Test: `nCat=1`, `trplMomIce=T`, `liqFrac=T` *(initial — pre-bugfix)*
  - Large apparent errors were recorded; root cause identified in 2.12

- [x] **2.10** Test: `nCat=2`, `trplMomIce=T`, `liqFrac=T` *(initial — pre-bugfix)*

- [x] **2.11** Full configuration sweep — all nCat × liqFrac combinations for 3-moment *(pre-bugfix)*
  - Ran via `full_testing/run_pp_tests`; new config files added: `1TF.txt`, `2TF.txt`, `3TF.txt`

- [x] **2.12** **Bug fix** — missing `Zitot_adv → Zitot_mp` conversion in `_get_diag_3mom()`
  - Root cause: `out_p3.dat` stores `Zitot` in the *advected* form `sqrt(Zi·Ni)`, but
    the LUT lookup requires the *microphysics* form `Zi²/Ni`.
    The conversion `zitot = zitot**2 / nitot` (Morrison et al. 2016, MWR) — applied at
    the top of `p3_main` — was missing from the post-processor.
  - Fix: added conversion in `_get_diag_3mom()` before the bounds clamp and LUT lookup.
  - Post-fix verification (all configs, `full_testing/run_pp_tests`):

  | Config      | diag_rhoi median rel err | diag_rhoi max rel err | diag_di median rel err | diag_di max rel err |
  |-------------|--------------------------|----------------------|------------------------|---------------------|
  | 1FF (2-mom) | 0.029%                   | 6.5%                 | 0.024%                 | 5.5%                |
  | 1TF         | 0.026%                   | 8.5%                 | 0.016%                 | 7.8%                |
  | 1TT         | 0.026%                   | 4.3%                 | 0.014%                 | 21.6%               |
  | 2TF cat 1   | 0.023%                   | 3.1%                 | 0.013%                 | 22.5%               |
  | 2TF cat 2   | 0.007%                   | 4.9%                 | 0.010%                 | 24.0%               |
  | 2TT cat 1   | 0.021%                   | 3.1%                 | 0.013%                 | 2.3%                |
  | 2TT cat 2   | 0.007%                   | 7.8%                 | 0.009%                 | 24.0%               |
  | 3TF cat 1   | 0.024%                   | 8.3%                 | 0.013%                 | 41.1%               |
  | 3TF cat 2   | 0.018%                   | 4.6%                 | 0.013%                 | 21.4%               |
  | 3TF cat 3   | 0.017%                   | 6.8%                 | 0.012%                 | 6.1%                |
  | 3TT cat 1   | 0.023%                   | 3.1%                 | 0.013%                 | 2.0%                |
  | 3TT cat 2   | 0.016%                   | 4.6%                 | 0.012%                 | 6.5%                |
  | 3TT cat 3   | 0.018%                   | 2.4%                 | 0.012%                 | 3.9%                |

  Median relative errors are ~0.01–0.03% for all configurations — consistent with the
  2-moment ASCII precision baseline. Max errors are isolated outliers at levels where
  `apply_mui_bounds_to_zi` changed `Zitot` (see note below).

**Note — residual max errors in 3-moment mode:**

After the `Zitot_adv → Zitot_mp` conversion fix, a small number of grid points still
show elevated max relative errors (up to ~41% for `diag_di` in a few configs). These
are the levels where `apply_mui_bounds_to_zi` in `p3_main` actively clamped `Zitot_mp`
to keep `mu_i` in bounds. The model's stored diagnostics were computed with the
pre-clamp `Zitot_mp`, while the post-processor reconstructs `Zitot_mp` from the
post-clamp `Zitot_adv` in `out_p3.dat`. At the vast majority of levels the bounds are
not active, so the median error is at the ASCII precision floor.

This is a file-format limitation, not a physics error. It can be resolved by storing
pre-clamp `Zitot_adv` (or the pre-clamp `Zitot_mp`) as an additional output column in
`cld1d` — deferred to a later phase if needed.

---

### Phase 3 — Input generalization

- [x] **3.1** Python API — accept NumPy arrays directly (no file I/O)
  - Decouple `P3LUT.get_diag()` from file reading (already mostly done)
  - Added `compute(arrays, cfg, lut_dir, lut)` top-level function
  - `main()` refactored to build arrays dict and call `compute()`
  - Test script: `postproc/test_compute.py`

- [ ] **3.2** netCDF input — read prognostic fields from `.nc` file
  - Use `netCDF4` or `xarray`
  - Auto-detect variable names (or accept a name mapping dict)
  - Handle units/scaling if needed

- [ ] **3.3** GEM FST/RPN input
  - Use `rpnpy` (CMC Python library) to read FST files
  - Map GEM variable names to P3 prognostic variable names
  - Handle grid and level structure

---

### Phase 4 — Additional diagnostics

- [ ] **4.1** `diag_effi` — effective ice radius (LUT index 6 for 2-mom, same for 3-mom)
- [ ] **4.2** `diag_vmi`  — mass-weighted fall speed (LUT index 2, needs `rhofaci`)
- [ ] **4.3** `diag_ZET`  — equivalent reflectivity (LUT index 9/13, needs air density)
- [ ] **4.4** Others as identified

---

## How to resume a session

1. Checkout the working branch:
   ```bash
   git checkout dev-postprocessor
   cd kin1d
   ```
2. Check this file (`PROJECT_PLAN.md`) for the first unchecked item
3. Read `CONTEXT.md` for technical reference (LUT dimensions, column layouts, constants)
4. The script to modify is `postproc/p3_postproc.py`
5. Run verification after each phase:
   ```bash
   ./execld                               # generate out_p3.dat
   python3 postproc/p3_postproc.py        # run post-processor
   ```

---

## Notes

- All code changes stay on branch `dev-postprocessor`; merge to `main` only at milestones
- Accuracy target: median relative error < 0.1% vs model output (baseline ~0.03%)
- Differences vs model expected due to ASCII precision loss in `out_p3.dat` (e13.4, ~5 sig figs)
- Binary input (Phase 3) will yield much smaller differences
