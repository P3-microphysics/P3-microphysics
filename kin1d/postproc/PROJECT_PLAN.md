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

- [x] **2.9** Test: `nCat=1`, `trplMomIce=T`, `liqFrac=T`
  - See **Known limitation** below — large apparent errors are expected and explained

- [x] **2.10** Test: `nCat=2`, `trplMomIce=T`, `liqFrac=T`

**Known limitation — Zitot time-lag in 3-moment mode:**

The `out_p3.dat` stores `Zitot` *after* `apply_mui_bounds_to_zi` has been applied
inside `p3_main`. However, `diag_rhoi` is computed earlier in the same timestep
using the *pre-bounds* `Zitot`. The post-processor reads the post-bounds `Zitot`
and therefore computes `diag_rhoi`/`diag_di` consistent with what the **next
timestep** would use — not the diagnostic stored alongside in the same row.

Consequence: comparison statistics against model `diag_rhoi`/`diag_di` show large
apparent errors (~30–100%) for 3-moment runs. The post-processor physics are
correct; the mismatch is a timing/output artefact in the ASCII file.

Options to resolve (deferred):
- Write pre-bounds `Zitot` as an additional column in `cld1d` output (Option B)
- Accept limitation for offline use (current approach)

---

### Phase 3 — Input generalization

- [ ] **3.1** Python API — accept NumPy arrays directly (no file I/O)
  - Decouple `P3LUT.get_diag()` from file reading (already mostly done)
  - Add a `compute(arrays_dict, config)` top-level function for programmatic use

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
   git checkout dev-tools_refactor
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

- All code changes stay on branch `dev-tools_refactor`; merge to `main` only at milestones
- Accuracy target: median relative error < 0.1% vs model output (baseline ~0.03%)
- Differences vs model expected due to ASCII precision loss in `out_p3.dat` (e13.4, ~5 sig figs)
- Binary input (Phase 3) will yield much smaller differences
