# P3 Post-Processor — Technical Context

A quick-reference document for resuming work on `p3_postproc.py` in a new session.
See `PROJECT_PLAN.md` for the current phase status and next steps.

---

## Repository layout (branch: dev-tools_refactor)

```
kin1d/
  postproc/
    p3_postproc.py         ← main post-processor script
    PROJECT_PLAN.md        ← project roadmap / phase status
    CONTEXT.md             ← this file
  src/
    p3/microphy_p3.f90     ← P3 source (key functions listed below)
    kin1d/cld1d.f90        ← kin1d driver (output format defined here)
    LUTs/create_p3_lookupTable_*.f90
  lookup_tables/           ← decompressed LUT files (gitignored)
    p3_lookupTable_1.dat-v6.9-2momI
    p3_lookupTable_1.dat-v6.9-3momI
    p3_lookupTable_3.dat-v1.4
  p3_config.txt            ← runtime config (nCat, trplMomIce, liqFrac)
  Makefile                 ← compiles from src/ subdirs directly
```

---

## Fortran → Python function mapping

| Fortran (microphy_p3.f90)              | Python (p3_postproc.py)             | Status   |
|----------------------------------------|-------------------------------------|----------|
| `calc_bulkRhoRime()`                   | `P3LUT.calc_bulk_rho_rime()`        | Done     |
| `find_lookupTable_indices_1a()`        | `P3LUT._find_indices()`             | Done     |
| `proc_from_LUT_main2mom()`             | `P3LUT._interp_2mom()`              | Done     |
| `find_lookupTable_indices_3a()`        | `P3LUT._find_indices_3a()`          | Phase 2  |
| `find_lookupTable_indices_1c()`        | `P3LUT._find_indices_1c()`          | Phase 2  |
| `proc_from_LUT_3()` (itab_3mom_mui)    | `P3LUT._interp_3mom_mui()`          | Phase 2  |
| `proc_from_LUT_main3mom()` (itab_3mom) | `P3LUT._interp_3mom()`              | Phase 2  |
| `get_mui_rhoi()`                       | incorporated in `get_diag_3mom()`   | Phase 2  |

---

## LUT array dimensions

| Array           | Shape                     | Contents                        |
|-----------------|---------------------------|---------------------------------|
| `itab`          | [5, 4, 4, 50, 19]         | 2-mom main LUT                  |
| `itab_3mom`     | [11, 5, 4, 4, 50, 29]     | 3-mom main LUT                  |
| `itab_3mom_mui` | [80, 5, 4, 4, 50, 2]      | 3-mom mu_i / rho_i LUT          |

Dimension order (all): `[zz or zq, densize, rimsize, liqsize, isize, process]`  
(2-mom has no zz dimension)

---

## LUT file structures

**`p3_lookupTable_1.dat-v6.9-2momI`** (2-moment):
- Line 1: header; Line 2: blank; Lines 3–: data
- Loop order (outer→inner): `jj → ii → ll → i`
- Per `(jj,ii,ll)` block: 50 main rows + 50×30 collection rows = 1,550 rows/block
- Columns per main row: `4 dummy | 19 values | 2 dummy` = **25 cols**
- Total lines: **124,002**

**`p3_lookupTable_1.dat-v6.9-3momI`** (3-moment):
- Line 1: header; Line 2: blank; Lines 3–: data
- Loop order (outer→inner): `zz → jj → ii → ll → i`
- Per `(zz,jj,ii,ll)` block: 50 main rows + 50×30 collection rows = 1,550 rows/block
- Columns per main row: `5 dummy | 29 values | 2 dummy` = **36 cols**
- Total lines: **1,364,002**

**`p3_lookupTable_3.dat-v1.4`** (mu_i LUT):
- Line 1: header; Line 2: blank; Lines 3–: data
- Loop order (outer→inner): `zq → jj → ii → ll → i`
- No collection table rows
- Columns per row: `5 dummy | 2 values | 2 dummy` = **9 cols**
- Total lines: **320,002**

---

## LUT process indices (1-based, matching Fortran)

**2-moment `itab`** (tabsize=19):

| Index | Quantity              |
|-------|-----------------------|
| 11    | `diag_di`  [m]        |
| 12    | `diag_rhoi` [kg m⁻³]  |

**3-moment `itab_3mom`** (tabsize_3mom=29):

| Index | Quantity              |
|-------|-----------------------|
| 11    | `diag_di`  [m]        |
| —     | `diag_rhoi` comes from `itab_3mom_mui` index 2, NOT from `itab_3mom` |

**`itab_3mom_mui`** (2 values):

| Index | Quantity              |
|-------|-----------------------|
| 1     | `mu_i`                |
| 2     | `diag_rhoi` [kg m⁻³]  |

---

## out_p3.dat column layout

**Prefix** (columns 0–14, same for all nCat):

| Col | Variable      | Col | Variable       |
|-----|---------------|-----|----------------|
| 0   | z             | 8   | Nc             |
| 1   | w             | 9   | Nr             |
| 2   | prt_liq       | 10  | Σ Qitot        |
| 3   | prt_sol       | 11  | Σ Nitot        |
| 4   | diag_ZET      | 12  | rime fraction  |
| 5   | T − T₀        | 13  | liq fraction   |
| 6   | Qc            | 14  | Drm (rain)     |
| 7   | Qr            |     |                |

**Per-category block** (8 columns each, starting at col 15):

| Offset | Variable    |
|--------|-------------|
| 0      | Qitot       |
| 1      | Qirim       |
| 2      | Qiliq       |
| 3      | Nitot       |
| 4      | Birim       |
| 5      | Zitot       |
| 6      | diag_rhoi   |
| 7      | diag_di     |

Total columns: `nCat=1` → 23, `nCat=2` → 31, `nCat=3` → 39

**Post-processor output (`out_pp.dat`):** same as `out_p3.dat` + `2×nCat` appended columns:
`pp_rhoi_cat1, pp_di_cat1 [, pp_rhoi_cat2, pp_di_cat2, ...]`

---

## Key constants (from microphy_p3.f90)

```python
QSMALL       = 1.e-14   # min ice mass mixing ratio     [kg kg⁻¹]
NSMALL       = 1.e-16   # min ice number concentration  [kg⁻¹]
RHO_RIME_MIN = 50.      # min bulk rime density         [kg m⁻³]
RHO_RIME_MAX = 900.     # max bulk rime density         [kg m⁻³]

# 2-moment LUT
ISIZE    = 50;  RIMSIZE  = 4;  LIQSIZE = 4;  DENSIZE = 5;  TABSIZE = 19;  RCOLLSIZE = 30

# 3-moment additions
ZSIZE = 11;  ZQSIZE = 80;  TABSIZE_3MOM = 29
```

---

## Index formulas

```
# size index (both 2-mom and 3-mom):
dum1  = (log10(Qitot/Nitot) + 18.) * 3.444606 - 10.    → clamped to [1, ISIZE]

# rime fraction index:
dum4  = (Qirim / (Qitot - Qiliq)) * 3. + 1.            → clamped to [1, RIMSIZE]

# liquid fraction index (= 1.0 when liqFrac=F):
dum7  = (Qiliq / Qitot) * 3. + 1.                      → clamped to [1, LIQSIZE]

# density index (non-linear):
dum5  = (rho_p - 50.) * 0.005 + 1.   if rho_p <= 650
      = (rho_p - 650.) * 0.004 + 4.  if rho_p >  650  → clamped to [1, DENSIZE]

# 3-moment only:
dum8  = (log10(Zitot/Qitot) + 23.) * 3.10347652        → clamped to [1, ZQSIZE]  (Zi/Qi → dumzq)
dum6  = mu_i * 0.5 + 1.                                → clamped to [1, ZSIZE]   (mu_i → dumzz)
```

---

## Verification baseline

| Config                          | Quantity   | Median rel err | Max rel err |
|---------------------------------|------------|----------------|-------------|
| nCat=1, 2-mom, liqFrac=F (done) | diag_rhoi  | 0.029%         | 6.5%        |
| nCat=1, 2-mom, liqFrac=F (done) | diag_di    | 0.024%         | 5.5%        |

Max errors are due to ASCII precision loss in `out_p3.dat` (e13.4 ≈ 5 significant digits).
Binary input (Phase 3) will reduce these to near machine precision.
