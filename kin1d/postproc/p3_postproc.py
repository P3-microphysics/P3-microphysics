#!/usr/bin/env python3
"""
p3_postproc.py — P3 microphysics offline post-processor
=========================================================

Computes diagnostic ice quantities (bulk density, mean-mass diameter, etc.)
from P3 prognostic variables by replicating the LUT lookup performed inside
p3_main. Intended for offline use when only prognostic model output is
available (e.g., GEM users wanting diagnostics after the fact).

QUICK START (from ./kin1d/):

    python3 postproc/p3_postproc.py

This reads:
    - out_p3.dat           (model output from cld1d)
    - lookup_tables/       (unzipped LUT files)
    - p3_config.txt        (P3 configuration: nCat, trplMomIce, liqFrac)

And writes:
    - out_pp.dat           (comparison output file)

The output file has the same format as out_p3.dat with 6 extra columns per
ice category appended (block of 6 per cat, cat 1 first):
    col N+0: pp_rhoi  [kg m-3]  ← bulk ice density
    col N+1: pp_di    [m]       ← mean-mass diameter
    col N+2: pp_effi  [m]       ← effective radius
    col N+3: pp_vmi   [m s-1]   ← mass-weighted fall speed (rhofaci=1; surface only)
    col N+4: pp_lami  [m-1]     ← PSD slope parameter lambda_i
    col N+5: pp_mui   [-]       ← PSD shape parameter mu_i

The original model diag_rhoi and diag_di columns are preserved, allowing
direct side-by-side comparison. Expected accuracy: median relative error ~0.03%.

COMMAND-LINE OPTIONS:

    --input FILE           Input file (default: out_p3.dat)
    --lut-dir DIR          LUT directory (default: lookup_tables)
    --config FILE          P3 config file (default: p3_config.txt)
    --output FILE          Output file (default: out_pp.dat)

SUPPORTED CONFIGURATIONS (v0.4):

    nCat=1–3, trplMomIce=.false., liqFrac=.false.  (2-moment, no liquid)
    nCat=1–3, trplMomIce=.false., liqFrac=.true.   (2-moment, with liquid fraction)
    nCat=1–3, trplMomIce=.true.,  liqFrac=.false.  (3-moment, no liquid)
    nCat=1–3, trplMomIce=.true.,  liqFrac=.true.   (3-moment, with liquid fraction)

PYTHON API:

    import sys
    sys.path.insert(0, 'postproc')   # if running from kin1d/
    from p3_postproc import compute, P3Config

    # Minimal example — 2-moment, single category, no liquid fraction:
    import numpy as np
    arrays = {
        'Qitot': np.array([...]),   # [kg kg-1]  shape (n_levels,) or (n_levels, n_cat)
        'Qirim': np.array([...]),   # [kg kg-1]
        'Qiliq': np.array([...]),   # [kg kg-1]  (zeros if liqFrac=F)
        'Nitot': np.array([...]),   # [kg-1]
        'Birim': np.array([...]),   # [m3 kg-1]
        'Zitot': np.array([...]),   # [m6 kg-1]  (zeros if trplMomIce=F)
    }
    cfg = P3Config(n_cat=1, triple_moment=False, liq_frac=False)
    result = compute(arrays, cfg, lut_dir='lookup_tables')
    pp_rhoi = result['pp_rhoi']   # shape (n_levels, n_cat)  [kg m-3]
    pp_di   = result['pp_di']     # shape (n_levels, n_cat)  [m]

    # Or load config from file:
    cfg = P3Config.from_file('p3_config.txt')
    result = compute(arrays, cfg, lut_dir='lookup_tables')

    Notes
    -----
    - Arrays may be 1-D (n_levels,) for n_cat=1, or 2-D (n_levels, n_cat).
    - 'Zitot' is required in the dict but ignored when triple_moment=False.
    - 'Qiliq' is required in the dict but treated as zero when liq_frac=False.
    - A P3LUT instance may be reused across compute() calls (pass via lut=):
          lut = P3LUT(lut_dir='lookup_tables', cfg=cfg)
          result = compute(arrays, cfg, lut=lut)

DEPENDENCIES:

    Python 3.6+
    NumPy   (pip3 install numpy)

See README.md in this directory for more details.
"""

from __future__ import annotations
import argparse
import numpy as np
from dataclasses import dataclass
from pathlib import Path


# ---------------------------------------------------------------------------
# Physical / numerical constants  (verbatim from microphy_p3.f90)
# ---------------------------------------------------------------------------
QSMALL       = 1.e-14   # minimum ice mass mixing ratio      [kg kg-1]
NSMALL       = 1.e-16   # minimum ice number concentration   [kg-1]
RHO_RIME_MIN = 50.      # minimum bulk rime density          [kg m-3]
RHO_RIME_MAX = 900.     # maximum bulk rime density          [kg m-3]
ZSMALL       = 1.e-35   # minimum Zitot (reflectivity moment)
ZLARGE       = 1.       # maximum Zitot (reflectivity moment)

# 2-moment LUT dimensions  (integer, parameter :: ... in microphy_p3.f90)
ISIZE     = 50    # particle size bins
RIMSIZE   =  4    # rime mass fraction bins
LIQSIZE   =  4    # liquid mass fraction bins
DENSIZE   =  5    # bulk rime density bins
TABSIZE   = 19    # number of process quantities in 2-moment table
RCOLLSIZE = 30    # ice-rain collection bins (file structure only; not used here)

# 3-moment additions
ZSIZE        = 11    # mu_i bins (outer dim of itab_3mom, = dum6/dumzz)
ZQSIZE       = 80    # Zi/Qi bins (outer dim of itab_3mom_mui, = dum8/dumzq)
TABSIZE_3MOM = 29    # number of process quantities in 3-moment table

# Process indices in the 2-moment LUT (1-based, matching Fortran)
IDX_DIAG_VMI  =  2    # mass-weighted fall speed  [m s-1]  (× rhofaci; here rhofaci=1)
IDX_DIAG_EFFI =  6    # effective radius          [m]
IDX_DIAG_DI   = 11    # mean-mass diameter        [m]
IDX_DIAG_RHOI = 12    # bulk ice density          [kg m-3]
IDX_DIAG_LAMI = 13    # lambda_i  (PSD slope)     [m-1]
IDX_DIAG_MUI  = 14    # mu_i      (PSD shape)     [-]

# 3-moment mu_i LUT indices (1-based, itab_3mom_mui)
IDX_3MOM_MUI  = 1    # mu_i parameter
IDX_3MOM_RHOI = 2    # diag_rhoi [kg m-3]

# 3-moment main LUT indices (1-based, itab_3mom)
IDX_3MOM_VMI  =  2   # mass-weighted fall speed  [m s-1]  (× rhofaci; here rhofaci=1)
IDX_3MOM_EFFI =  6   # effective radius          [m]
IDX_3MOM_DI   = 11   # diag_di   [m]
IDX_3MOM_LAMI = 14   # lambda_i  [m-1]
# Note: mu_i (index 15 in itab_3mom) is not fetched here — it is taken directly
# from the IDX_3MOM_MUI intermediate already computed via itab_3mom_mui.


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass
class P3Config:
    """P3 runtime configuration, parsed from p3_config.txt."""
    n_cat:         int  = 1
    triple_moment: bool = False
    liq_frac:      bool = False

    @classmethod
    def from_file(cls, path: str) -> P3Config:
        """
        Parse a p3_config.txt file.

        Expected format (value before the ! comment):
            1        ! nCat
            .false.  ! trplMomIce
            .true.   ! liqFrac
        """
        cfg = cls()
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if '!' not in line:
                    continue
                val_str, _, key_str = line.partition('!')
                val = val_str.strip().lower()
                key = key_str.strip().lower()
                if key == 'ncat':
                    cfg.n_cat = int(val)
                elif key == 'trplmomice':
                    cfg.triple_moment = val == '.true.'
                elif key == 'liqfrac':
                    cfg.liq_frac = val == '.true.'
        return cfg


# ---------------------------------------------------------------------------
# LUT loader and interpolator
# ---------------------------------------------------------------------------
class P3LUT:
    """
    Loads a P3 lookup table and computes diagnostic ice quantities.

    For 2-moment ice, mirrors:
        calc_bulkRhoRime()             in microphy_p3.f90
        find_lookupTable_indices_1a()  "
        proc_from_LUT_main2mom()       "

    For 3-moment ice, additionally mirrors:
        find_lookupTable_indices_3a()  "
        find_lookupTable_indices_1c()  "
        proc_from_LUT_3()              "  (itab_3mom_mui)
        proc_from_LUT_main3mom()       "  (itab_3mom)
        get_mui_rhoi()                 "
    """

    def __init__(self, lut_dir: str, cfg: P3Config):
        self.cfg = cfg
        lut_dir  = Path(lut_dir)
        if not cfg.triple_moment:
            candidates = list(lut_dir.glob('p3_lookupTable_1.dat-*2momI*'))
            if not candidates:
                raise FileNotFoundError(
                    f'No 2-moment LUT file found in {lut_dir}. '
                    f'Make sure the .gz files are decompressed.')
            self._load_2mom_lut(candidates[0])
        else:
            self._load_3mom_luts(lut_dir)

    # ------------------------------------------------------------------
    # 2-moment LUT loader
    # ------------------------------------------------------------------
    def _load_2mom_lut(self, path: Path):
        """
        Read the 2-moment LUT into self.itab[densize, rimsize, liqsize, isize, tabsize].

        File structure mirrors the Fortran p3_init read loop:
          - 1 header line + 1 blank line (skipped)
          - Loop (outer→inner): jj → ii → ll → i
          - Per (jj, ii, ll): isize rows × 25 cols: 4 dummy | 19 table values | 2 dummy
          - Per (jj, ii, ll): isize×rcollsize rows  (ice-rain collection data, skipped)
        """
        itab = np.zeros((DENSIZE, RIMSIZE, LIQSIZE, ISIZE, TABSIZE), dtype=np.float32)
        with open(path) as fh:
            fh.readline()   # header line
            fh.readline()   # blank line after header
            for jj in range(DENSIZE):
                for ii in range(RIMSIZE):
                    for ll in range(LIQSIZE):
                        # main table: isize rows
                        for i in range(ISIZE):
                            vals = list(map(float, fh.readline().split()))
                            itab[jj, ii, ll, i, :] = vals[4 : 4 + TABSIZE]
                        # collection table: skip isize × rcollsize rows
                        for _ in range(ISIZE * RCOLLSIZE):
                            fh.readline()
        self.itab = itab
        print(f'  Loaded 2-moment LUT: {path.name}  shape={itab.shape}')

    # ------------------------------------------------------------------
    # 3-moment LUT loaders
    # ------------------------------------------------------------------
    def _load_3mom_luts(self, lut_dir: Path):
        """Load both 3-moment LUT files."""
        # main 3-mom LUT (itab_3mom)
        candidates = list(lut_dir.glob('p3_lookupTable_1.dat-*3momI*'))
        if not candidates:
            raise FileNotFoundError(
                f'No 3-moment LUT file found in {lut_dir}. '
                f'Make sure the .gz files are decompressed.')
        self._load_3mom_main_lut(candidates[0])

        # mu_i / rho_i LUT (itab_3mom_mui)
        candidates = list(lut_dir.glob('p3_lookupTable_3.dat-*'))
        if not candidates:
            raise FileNotFoundError(
                f'No p3_lookupTable_3 file found in {lut_dir}. '
                f'Make sure the .gz files are decompressed.')
        self._load_3mom_mui_lut(candidates[0])

    def _load_3mom_main_lut(self, path: Path):
        """
        Read the 3-moment main LUT into
        self.itab_3mom[zsize, densize, rimsize, liqsize, isize, tabsize_3mom].

        File structure mirrors the Fortran p3_init read loop:
          - 1 header line + 1 blank line (skipped)
          - Loop (outer→inner): zz → jj → ii → ll → i
          - Per (zz, jj, ii, ll): isize rows × 36 cols:
              5 dummy | 29 table values | 2 dummy
          - Per (zz, jj, ii, ll): isize×rcollsize collection rows (4 dummy + 3 vals, skipped)
        """
        shape = (ZSIZE, DENSIZE, RIMSIZE, LIQSIZE, ISIZE, TABSIZE_3MOM)
        itab = np.zeros(shape, dtype=np.float32)
        with open(path) as fh:
            fh.readline()   # header line
            fh.readline()   # blank line after header
            for zz in range(ZSIZE):
                for jj in range(DENSIZE):
                    for ii in range(RIMSIZE):
                        for ll in range(LIQSIZE):
                            # main table: isize rows, 5 dummy + 29 values + 2 dummy
                            for i in range(ISIZE):
                                vals = list(map(float, fh.readline().split()))
                                itab[zz, jj, ii, ll, i, :] = vals[5 : 5 + TABSIZE_3MOM]
                            # collection table: isize × rcollsize rows (skip)
                            for _ in range(ISIZE * RCOLLSIZE):
                                fh.readline()
        self.itab_3mom = itab
        print(f'  Loaded 3-moment main LUT: {path.name}  shape={itab.shape}')

    def _load_3mom_mui_lut(self, path: Path):
        """
        Read the mu_i/rho_i LUT into
        self.itab_3mom_mui[zqsize, densize, rimsize, liqsize, isize, 2].

        File structure mirrors the Fortran p3_init read loop:
          - 1 header line + 1 blank line (skipped)
          - Loop (outer→inner): zq → jj → ii → ll → i
          - Per row: 9 cols: 5 dummy | 2 values | 2 dummy
          - No collection table rows
        """
        shape = (ZQSIZE, DENSIZE, RIMSIZE, LIQSIZE, ISIZE, 2)
        itab = np.zeros(shape, dtype=np.float32)
        with open(path) as fh:
            fh.readline()   # header line
            fh.readline()   # blank line after header
            for zq in range(ZQSIZE):
                for jj in range(DENSIZE):
                    for ii in range(RIMSIZE):
                        for ll in range(LIQSIZE):
                            for i in range(ISIZE):
                                vals = list(map(float, fh.readline().split()))
                                itab[zq, jj, ii, ll, i, :] = vals[5:7]
        self.itab_3mom_mui = itab
        print(f'  Loaded 3-moment mu_i LUT: {path.name}  shape={itab.shape}')

    # ------------------------------------------------------------------
    # Bulk rime density
    # ------------------------------------------------------------------
    @staticmethod
    def calc_bulk_rho_rime(
            qi_tot: float, qi_rim: float, qi_liq: float, bi_rim: float
    ) -> tuple[float, float, float]:
        """
        Compute bulk rime density and apply consistency limits.

        Verbatim translation of calc_bulkRhoRime() from microphy_p3.f90.

        Returns
        -------
        rho_rime : float  [kg m-3]
        qi_rim   : float  [kg kg-1]  (possibly adjusted)
        bi_rim   : float  [m3 kg-1]  (possibly adjusted)
        """
        if bi_rim >= 1.e-15:
            rho_rime = qi_rim / bi_rim
            if rho_rime < RHO_RIME_MIN:
                rho_rime = RHO_RIME_MIN
                bi_rim   = qi_rim / rho_rime
            elif rho_rime > RHO_RIME_MAX:
                rho_rime = RHO_RIME_MAX
                bi_rim   = qi_rim / rho_rime
        else:
            qi_rim = 0.; bi_rim = 0.; rho_rime = 0.

        if qi_rim < QSMALL:
            qi_rim = 0.; bi_rim = 0.
        elif qi_rim > (qi_tot - qi_liq) and rho_rime > 0.:
            qi_rim = qi_tot - qi_liq
            bi_rim = qi_rim / rho_rime

        return rho_rime, qi_rim, bi_rim

    # ------------------------------------------------------------------
    # Index finding — shared (2-mom and 3-mom)
    # ------------------------------------------------------------------
    @staticmethod
    def _find_indices(
            qi_tot: float, ni_tot: float, qi_rim: float,
            qi_liq: float, rho_p: float
    ) -> tuple[int, int, int, int, float, float, float, float]:
        """
        Compute LUT indices from prognostic ice variables.

        Translation of find_lookupTable_indices_1a() from microphy_p3.f90.
        All integer indices are 1-based (Fortran convention) and clamped to
        valid interpolation range.

        Returns
        -------
        dumi, dumjj, dumii, dumll : int    lower-bound integer indices (1-based)
        dum1, dum5, dum4, dum7   : float  continuous indices (clamped to bounds)
        """
        # size index  (from qi_tot / ni_tot ratio)
        dum1  = (np.log10(qi_tot / ni_tot) + 18.) * 3.444606 - 10.
        dum1  = float(np.clip(dum1, 1., ISIZE))
        dumi  = int(np.clip(int(dum1), 1, ISIZE - 1))

        # rime fraction index
        qi_dry = qi_tot - qi_liq
        frim   = (qi_rim / qi_dry) if qi_dry > 0. else 0.
        dum4   = frim * 3. + 1.
        dum4   = float(np.clip(dum4, 1., RIMSIZE))
        dumii  = int(np.clip(int(dum4), 1, RIMSIZE - 1))

        # liquid fraction index
        dum7  = (qi_liq / qi_tot) * 3. + 1.
        dum7  = float(np.clip(dum7, 1., LIQSIZE))
        dumll = int(np.clip(int(dum7), 1, LIQSIZE - 1))

        # density index  (non-linear spacing, matches Fortran exactly)
        if rho_p <= 650.:
            dum5 = (rho_p - 50.)  * 0.005 + 1.
        else:
            dum5 = (rho_p - 650.) * 0.004 + 4.
        dum5  = float(np.clip(dum5, 1., DENSIZE))
        dumjj = int(np.clip(int(dum5), 1, DENSIZE - 1))

        return dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7

    # ------------------------------------------------------------------
    # Index finding — 3-moment only
    # ------------------------------------------------------------------
    @staticmethod
    def _find_indices_3a(
            zi_tot: float, qi_tot: float
    ) -> tuple[int, float]:
        """
        Compute Zi/Qi index for the 3-moment mu_i LUT.

        Translation of find_lookupTable_indices_3a() from microphy_p3.f90.

        Returns
        -------
        dumzq : int    lower-bound integer index (1-based)
        dum8  : float  continuous index (clamped to [1, ZQSIZE])
        """
        dum8  = (np.log10(zi_tot / qi_tot) + 23.) * 3.10347652
        dumzq = int(dum8)
        dum8  = float(np.clip(dum8, 1., ZQSIZE))
        dumzq = int(np.clip(dumzq, 1, ZQSIZE - 1))
        return dumzq, dum8

    @staticmethod
    def _find_indices_1c(
            mu_i: float
    ) -> tuple[int, float]:
        """
        Compute mu_i index for the 3-moment main LUT.

        Translation of find_lookupTable_indices_1c() from microphy_p3.f90.

        Returns
        -------
        dumzz : int    lower-bound integer index (1-based)
        dum6  : float  continuous index (clamped to [1, ZSIZE])
        """
        dum6  = mu_i * 0.5 + 1.
        dumzz = int(dum6)
        dum6  = float(np.clip(dum6, 1., ZSIZE))
        dumzz = int(np.clip(dumzz, 1, ZSIZE - 1))
        return dumzz, dum6

    # ------------------------------------------------------------------
    # 2-moment interpolation
    # ------------------------------------------------------------------
    def _interp_2mom(
            self,
            ind_f: int,
            dumi: int, dumjj: int, dumii: int, dumll: int,
            dum1: float, dum5: float, dum4: float, dum7: float
    ) -> float:
        """
        4D (or 3D) linear interpolation in itab.

        Verbatim translation of proc_from_LUT_main2mom() from microphy_p3.f90.

        Parameters
        ----------
        ind_f : 1-based process index (11 = diag_di, 12 = diag_rhoi)
        dumi, dumjj, dumii, dumll : 1-based lower-bound integer indices
        dum1  : continuous size index
        dum5  : continuous density index
        dum4  : continuous rime-fraction index
        dum7  : continuous liquid-fraction index
        """
        ind = ind_f - 1    # convert to 0-based for numpy
        t   = self.itab    # shape: [jj, ii, ll, i, ind]  (all 0-based)

        # convert Fortran 1-based lower-bound indices → 0-based
        jj = dumjj - 1;  jj1 = jj + 1
        ii = dumii - 1;  ii1 = ii + 1
        ll = dumll - 1;  ll1 = ll + 1
        i  = dumi  - 1;  i1  = i  + 1

        # interpolation fractions (continuous index minus its 1-based floor)
        fi  = dum1 - dumi    # size
        fii = dum4 - dumii   # rime fraction
        fjj = dum5 - dumjj   # density
        fll = dum7 - dumll   # liquid fraction

        if dum7 == 1. and dumll == 1:
            # ---- 3D path: liquid fraction dimension skipped ----
            # at density jj: interp in size, then rime
            iproc1 = t[jj,  ii,  ll, i,  ind] + fi * (t[jj,  ii,  ll, i1, ind] - t[jj,  ii,  ll, i,  ind])
            iproc2 = t[jj,  ii1, ll, i,  ind] + fi * (t[jj,  ii1, ll, i1, ind] - t[jj,  ii1, ll, i,  ind])
            tmp1   = iproc1 + fii * (iproc2 - iproc1)

            # at density jj+1: interp in size, then rime
            iproc1 = t[jj1, ii,  ll, i,  ind] + fi * (t[jj1, ii,  ll, i1, ind] - t[jj1, ii,  ll, i,  ind])
            iproc2 = t[jj1, ii1, ll, i,  ind] + fi * (t[jj1, ii1, ll, i1, ind] - t[jj1, ii1, ll, i,  ind])
            tmp2   = iproc1 + fii * (iproc2 - iproc1)

            return float(tmp1 + fjj * (tmp2 - tmp1))

        else:
            # ---- 4D path: liquid fraction dimension included ----
            # at density jj, rime ii:
            iproc1 = t[jj,  ii,  ll,  i, ind] + fi * (t[jj,  ii,  ll,  i1, ind] - t[jj,  ii,  ll,  i, ind])
            iproc2 = t[jj,  ii,  ll1, i, ind] + fi * (t[jj,  ii,  ll1, i1, ind] - t[jj,  ii,  ll1, i, ind])
            tmp1   = iproc1 + fll * (iproc2 - iproc1)

            # at density jj, rime ii+1:
            iproc3 = t[jj,  ii1, ll,  i, ind] + fi * (t[jj,  ii1, ll,  i1, ind] - t[jj,  ii1, ll,  i, ind])
            iproc4 = t[jj,  ii1, ll1, i, ind] + fi * (t[jj,  ii1, ll1, i1, ind] - t[jj,  ii1, ll1, i, ind])
            tmp2   = iproc3 + fll * (iproc4 - iproc3)
            tmp3   = tmp1 + fii * (tmp2 - tmp1)   # interp in rime  → result at jj

            # at density jj+1, rime ii:
            iproc1 = t[jj1, ii,  ll,  i, ind] + fi * (t[jj1, ii,  ll,  i1, ind] - t[jj1, ii,  ll,  i, ind])
            iproc2 = t[jj1, ii,  ll1, i, ind] + fi * (t[jj1, ii,  ll1, i1, ind] - t[jj1, ii,  ll1, i, ind])
            tmp1   = iproc1 + fll * (iproc2 - iproc1)

            # at density jj+1, rime ii+1:
            iproc3 = t[jj1, ii1, ll,  i, ind] + fi * (t[jj1, ii1, ll,  i1, ind] - t[jj1, ii1, ll,  i, ind])
            iproc4 = t[jj1, ii1, ll1, i, ind] + fi * (t[jj1, ii1, ll1, i1, ind] - t[jj1, ii1, ll1, i, ind])
            tmp2   = iproc3 + fll * (iproc4 - iproc3)
            tmp4   = tmp1 + fii * (tmp2 - tmp1)   # interp in rime  → result at jj+1

            return float(tmp3 + fjj * (tmp4 - tmp3))   # interp in density

    # ------------------------------------------------------------------
    # 3-moment interpolation — mu_i / rho_i LUT  (itab_3mom_mui)
    # ------------------------------------------------------------------
    def _interp_3mom_mui(
            self,
            ind_f: int,
            dumzq: int, dumjj: int, dumii: int, dumll: int, dumi: int,
            dum1: float, dum4: float, dum5: float, dum7: float, dum8: float
    ) -> float:
        """
        5D linear interpolation in itab_3mom_mui.

        Verbatim translation of proc_from_LUT_3() from microphy_p3.f90.

        Parameters
        ----------
        ind_f  : 1-based process index (1 = mu_i, 2 = diag_rhoi)
        dumzq  : 1-based lower-bound Zi/Qi index
        dumjj, dumii, dumll, dumi : 1-based lower-bound density/rime/liq/size indices
        dum1   : continuous size index
        dum4   : continuous rime-fraction index
        dum5   : continuous density index
        dum7   : continuous liquid-fraction index
        dum8   : continuous Zi/Qi index
        """
        ind = ind_f - 1
        t   = self.itab_3mom_mui   # shape: [zq, jj, ii, ll, i, ind]  (0-based)

        # convert to 0-based
        zq = dumzq - 1;  zq1 = zq + 1
        jj = dumjj - 1;  jj1 = jj + 1
        ii = dumii - 1;  ii1 = ii + 1
        ll = dumll - 1;  ll1 = ll + 1
        i  = dumi  - 1;  i1  = i  + 1

        fi  = dum1 - dumi
        fii = dum4 - dumii
        fjj = dum5 - dumjj
        fll = dum7 - dumll
        fzq = dum8 - dumzq

        def _interp_4d(zq_):
            """4D or 3D interpolation at a fixed zq index."""
            if dum7 == 1. and dumll == 1:
                # 3D path: skip liquid fraction
                dp1  = t[zq_, jj,  ii,  ll, i,  ind] + fi * (t[zq_, jj,  ii,  ll, i1, ind] - t[zq_, jj,  ii,  ll, i,  ind])
                dp2  = t[zq_, jj,  ii1, ll, i,  ind] + fi * (t[zq_, jj,  ii1, ll, i1, ind] - t[zq_, jj,  ii1, ll, i,  ind])
                gp1  = dp1 + fii * (dp2 - dp1)
                dp1  = t[zq_, jj1, ii,  ll, i,  ind] + fi * (t[zq_, jj1, ii,  ll, i1, ind] - t[zq_, jj1, ii,  ll, i,  ind])
                dp2  = t[zq_, jj1, ii1, ll, i,  ind] + fi * (t[zq_, jj1, ii1, ll, i1, ind] - t[zq_, jj1, ii1, ll, i,  ind])
                gp2  = dp1 + fii * (dp2 - dp1)
                return gp1 + fjj * (gp2 - gp1)
            else:
                # 4D path: include liquid fraction
                # jj, ii
                dp1  = t[zq_, jj,  ii,  ll,  i, ind] + fi * (t[zq_, jj,  ii,  ll,  i1, ind] - t[zq_, jj,  ii,  ll,  i, ind])
                dp2  = t[zq_, jj,  ii,  ll1, i, ind] + fi * (t[zq_, jj,  ii,  ll1, i1, ind] - t[zq_, jj,  ii,  ll1, i, ind])
                ip1  = dp1 + fll * (dp2 - dp1)
                # jj, ii+1
                dp1  = t[zq_, jj,  ii1, ll,  i, ind] + fi * (t[zq_, jj,  ii1, ll,  i1, ind] - t[zq_, jj,  ii1, ll,  i, ind])
                dp2  = t[zq_, jj,  ii1, ll1, i, ind] + fi * (t[zq_, jj,  ii1, ll1, i1, ind] - t[zq_, jj,  ii1, ll1, i, ind])
                ip2  = dp1 + fll * (dp2 - dp1)
                gp1  = ip1 + fii * (ip2 - ip1)
                # jj+1, ii
                dp1  = t[zq_, jj1, ii,  ll,  i, ind] + fi * (t[zq_, jj1, ii,  ll,  i1, ind] - t[zq_, jj1, ii,  ll,  i, ind])
                dp2  = t[zq_, jj1, ii,  ll1, i, ind] + fi * (t[zq_, jj1, ii,  ll1, i1, ind] - t[zq_, jj1, ii,  ll1, i, ind])
                ip1  = dp1 + fll * (dp2 - dp1)
                # jj+1, ii+1
                dp1  = t[zq_, jj1, ii1, ll,  i, ind] + fi * (t[zq_, jj1, ii1, ll,  i1, ind] - t[zq_, jj1, ii1, ll,  i, ind])
                dp2  = t[zq_, jj1, ii1, ll1, i, ind] + fi * (t[zq_, jj1, ii1, ll1, i1, ind] - t[zq_, jj1, ii1, ll1, i, ind])
                ip2  = dp1 + fll * (dp2 - dp1)
                gp2  = ip1 + fii * (ip2 - ip1)
                return gp1 + fjj * (gp2 - gp1)

        rproc1 = _interp_4d(zq)
        rproc2 = _interp_4d(zq1)
        return float(rproc1 + fzq * (rproc2 - rproc1))

    # ------------------------------------------------------------------
    # 3-moment interpolation — main LUT  (itab_3mom)
    # ------------------------------------------------------------------
    def _interp_3mom(
            self,
            ind_f: int,
            dumzz: int, dumjj: int, dumii: int, dumll: int, dumi: int,
            dum1: float, dum4: float, dum5: float, dum6: float, dum7: float
    ) -> float:
        """
        6D (or 5D) linear interpolation in itab_3mom.

        Verbatim translation of proc_from_LUT_main3mom() from microphy_p3.f90.

        Parameters
        ----------
        ind_f  : 1-based process index (11 = diag_di)
        dumzz  : 1-based lower-bound mu_i index
        dumjj, dumii, dumll, dumi : 1-based lower-bound density/rime/liq/size indices
        dum1   : continuous size index
        dum4   : continuous rime-fraction index
        dum5   : continuous density index
        dum6   : continuous mu_i index
        dum7   : continuous liquid-fraction index
        """
        ind = ind_f - 1
        t   = self.itab_3mom   # shape: [zz, jj, ii, ll, i, ind]  (0-based)

        # convert to 0-based
        zz = dumzz - 1;  zz1 = zz + 1
        jj = dumjj - 1;  jj1 = jj + 1
        ii = dumii - 1;  ii1 = ii + 1
        ll = dumll - 1;  ll1 = ll + 1
        i  = dumi  - 1;  i1  = i  + 1

        fi  = dum1 - dumi
        fii = dum4 - dumii
        fjj = dum5 - dumjj
        fll = dum7 - dumll
        fzz = dum6 - dumzz

        def _interp_4d(zz_):
            """4D or 3D interpolation at a fixed zz index."""
            if dum7 == 1. and dumll == 1:
                # 3D path: skip liquid fraction
                dp1  = t[zz_, jj,  ii,  ll, i,  ind] + fi * (t[zz_, jj,  ii,  ll, i1, ind] - t[zz_, jj,  ii,  ll, i,  ind])
                dp2  = t[zz_, jj,  ii1, ll, i,  ind] + fi * (t[zz_, jj,  ii1, ll, i1, ind] - t[zz_, jj,  ii1, ll, i,  ind])
                gp1  = dp1 + fii * (dp2 - dp1)
                dp1  = t[zz_, jj1, ii,  ll, i,  ind] + fi * (t[zz_, jj1, ii,  ll, i1, ind] - t[zz_, jj1, ii,  ll, i,  ind])
                dp2  = t[zz_, jj1, ii1, ll, i,  ind] + fi * (t[zz_, jj1, ii1, ll, i1, ind] - t[zz_, jj1, ii1, ll, i,  ind])
                gp2  = dp1 + fii * (dp2 - dp1)
                return gp1 + fjj * (gp2 - gp1)
            else:
                # 4D path: include liquid fraction
                # jj, ii
                dp1  = t[zz_, jj,  ii,  ll,  i, ind] + fi * (t[zz_, jj,  ii,  ll,  i1, ind] - t[zz_, jj,  ii,  ll,  i, ind])
                dp2  = t[zz_, jj,  ii,  ll1, i, ind] + fi * (t[zz_, jj,  ii,  ll1, i1, ind] - t[zz_, jj,  ii,  ll1, i, ind])
                ip1  = dp1 + fll * (dp2 - dp1)
                # jj, ii+1
                dp1  = t[zz_, jj,  ii1, ll,  i, ind] + fi * (t[zz_, jj,  ii1, ll,  i1, ind] - t[zz_, jj,  ii1, ll,  i, ind])
                dp2  = t[zz_, jj,  ii1, ll1, i, ind] + fi * (t[zz_, jj,  ii1, ll1, i1, ind] - t[zz_, jj,  ii1, ll1, i, ind])
                ip2  = dp1 + fll * (dp2 - dp1)
                gp1  = ip1 + fii * (ip2 - ip1)
                # jj+1, ii
                dp1  = t[zz_, jj1, ii,  ll,  i, ind] + fi * (t[zz_, jj1, ii,  ll,  i1, ind] - t[zz_, jj1, ii,  ll,  i, ind])
                dp2  = t[zz_, jj1, ii,  ll1, i, ind] + fi * (t[zz_, jj1, ii,  ll1, i1, ind] - t[zz_, jj1, ii,  ll1, i, ind])
                ip1  = dp1 + fll * (dp2 - dp1)
                # jj+1, ii+1
                dp1  = t[zz_, jj1, ii1, ll,  i, ind] + fi * (t[zz_, jj1, ii1, ll,  i1, ind] - t[zz_, jj1, ii1, ll,  i, ind])
                dp2  = t[zz_, jj1, ii1, ll1, i, ind] + fi * (t[zz_, jj1, ii1, ll1, i1, ind] - t[zz_, jj1, ii1, ll1, i, ind])
                ip2  = dp1 + fll * (dp2 - dp1)
                gp2  = ip1 + fii * (ip2 - ip1)
                return gp1 + fjj * (gp2 - gp1)

        rproc1 = _interp_4d(zz)
        rproc2 = _interp_4d(zz1)
        return float(rproc1 + fzz * (rproc2 - rproc1))

    # ------------------------------------------------------------------
    # Top-level diagnostic calls
    # ------------------------------------------------------------------
    def get_diag(
            self,
            qi_tot: float, qi_rim: float, qi_liq: float,
            ni_tot: float, bi_rim: float,
            zi_tot: float = 0.
    ) -> tuple[float, float, float, float, float, float]:
        """
        Compute ice diagnostics for a single grid point.

        Routes to the 2-moment or 3-moment path depending on cfg.triple_moment.
        Returns all-zeros when qi_tot < QSMALL, matching model behaviour.

        zi_tot (3-moment only): the *advected* form of the sixth moment,
            i.e. sqrt(Zi * Ni), as written to out_p3.dat by cld1d.
            The conversion to the microphysics form (Zi²/Ni) is applied
            internally in _get_diag_3mom() before the LUT lookup.

        Returns
        -------
        diag_rhoi, diag_di, diag_effi, diag_vmi, diag_lami, diag_mui

        Note: diag_vmi uses rhofaci = 1.0 (valid at the surface only).
        """
        if qi_tot < QSMALL:
            return 0., 0., 0., 0., 0., 0.

        if self.cfg.triple_moment:
            return self._get_diag_3mom(qi_tot, qi_rim, qi_liq, ni_tot, bi_rim, zi_tot)
        else:
            return self._get_diag_2mom(qi_tot, qi_rim, qi_liq, ni_tot, bi_rim)

    def _get_diag_2mom(
            self,
            qi_tot: float, qi_rim: float, qi_liq: float,
            ni_tot: float, bi_rim: float
    ) -> tuple[float, float, float, float, float, float]:
        """
        2-moment path: compute all ice diagnostics.

        Replicates the qi_not_small block in p3_main (microphy_p3.f90).

        Returns: diag_rhoi, diag_di, diag_effi, diag_vmi, diag_lami, diag_mui
        Note: diag_vmi uses rhofaci = 1.0 (valid at the surface only).
        """
        ni_tot = max(ni_tot, NSMALL)

        rho_p, qi_rim, bi_rim = self.calc_bulk_rho_rime(
            qi_tot, qi_rim, qi_liq, bi_rim)

        dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7 = self._find_indices(
            qi_tot, ni_tot, qi_rim, qi_liq, rho_p)

        args = (dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7)
        diag_rhoi = self._interp_2mom(IDX_DIAG_RHOI, *args)
        diag_di   = self._interp_2mom(IDX_DIAG_DI,   *args)
        diag_effi = self._interp_2mom(IDX_DIAG_EFFI, *args)
        diag_vmi  = self._interp_2mom(IDX_DIAG_VMI,  *args)   # rhofaci = 1.0
        diag_lami = self._interp_2mom(IDX_DIAG_LAMI, *args)
        diag_mui  = self._interp_2mom(IDX_DIAG_MUI,  *args)

        return diag_rhoi, diag_di, diag_effi, diag_vmi, diag_lami, diag_mui

    def _get_diag_3mom(
            self,
            qi_tot: float, qi_rim: float, qi_liq: float,
            ni_tot: float, bi_rim: float,
            zi_tot: float
    ) -> tuple[float, float, float, float, float, float]:
        """
        3-moment path: compute all ice diagnostics.

        zi_tot must be the *advected* form sqrt(Zi*Ni) as stored in out_p3.dat.
        It is converted to the microphysics form Zi²/Ni here before LUT lookup,
        mirroring the conversion at the top of p3_main (Morrison et al. 2016, MWR).

        Mirrors the get_mui_rhoi() call sequence in p3_main (microphy_p3.f90):
          1. calc_bulkRhoRime → rho_p
          2. find_lookupTable_indices_1a → shared indices
          3. Convert Zitot_adv → Zitot_mp; clamp; find_lookupTable_indices_3a → dumzq, dum8
          4. proc_from_LUT_3(1) → mu_i
          5. proc_from_LUT_3(2) → diag_rhoi
          6. find_lookupTable_indices_1c → dumzz, dum6
          7. proc_from_LUT_main3mom(2,6,11,14) → diag_vmi, diag_effi, diag_di, diag_lami

        Returns: diag_rhoi, diag_di, diag_effi, diag_vmi, diag_lami, diag_mui
        Note: diag_vmi uses rhofaci = 1.0 (valid at the surface only).
              diag_mui is taken directly from the mu_i intermediate (step 4).
        """
        ni_tot = max(ni_tot, NSMALL)

        rho_p, qi_rim, bi_rim = self.calc_bulk_rho_rime(
            qi_tot, qi_rim, qi_liq, bi_rim)

        dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7 = self._find_indices(
            qi_tot, ni_tot, qi_rim, qi_liq, rho_p)

        # Convert Zitot from the advected form sqrt(Zi*Ni) stored in out_p3.dat
        # to the microphysics form Zi²/Ni used by the LUT.
        # Mirrors the conversion at the top of p3_main (microphy_p3.f90 line ~2408):
        #   where (nitot>0.) zitot = zitot**2/nitot   ! Morrison et al. (2016), MWR
        # ni_tot >= NSMALL at this point, so no zero-division risk.
        zi_tot = zi_tot ** 2 / ni_tot

        # clamp Zitot_mp (mirrors 'Zi = max(zsmall, min(zlarge, Zi))' in get_mui_rhoi)
        zi_tot = max(ZSMALL, min(ZLARGE, zi_tot))

        # find Zi/Qi index and interpolate mu_i and rho_i from itab_3mom_mui
        dumzq, dum8 = self._find_indices_3a(zi_tot, qi_tot)

        mui_args = (dumzq, dumjj, dumii, dumll, dumi, dum1, dum4, dum5, dum7, dum8)
        mu_i      = self._interp_3mom_mui(IDX_3MOM_MUI,  *mui_args)
        diag_rhoi = self._interp_3mom_mui(IDX_3MOM_RHOI, *mui_args)

        # find mu_i index and interpolate remaining fields from itab_3mom
        dumzz, dum6 = self._find_indices_1c(mu_i)

        main_args = (dumzz, dumjj, dumii, dumll, dumi, dum1, dum4, dum5, dum6, dum7)
        diag_di   = self._interp_3mom(IDX_3MOM_DI,   *main_args)
        diag_effi = self._interp_3mom(IDX_3MOM_EFFI, *main_args)
        diag_vmi  = self._interp_3mom(IDX_3MOM_VMI,  *main_args)   # rhofaci = 1.0
        diag_lami = self._interp_3mom(IDX_3MOM_LAMI, *main_args)
        diag_mui  = mu_i   # already computed from itab_3mom_mui

        return diag_rhoi, diag_di, diag_effi, diag_vmi, diag_lami, diag_mui


# ---------------------------------------------------------------------------
# Column layout for out_p3.dat / out_pp.dat
# ---------------------------------------------------------------------------
# The output file from cld1d has a fixed prefix of 15 columns (cols 0-14)
# followed by one block of 8 columns per ice category:
#   col 0: Qitot   col 4: Birim
#   col 1: Qirim   col 5: Zitot
#   col 2: Qiliq   col 6: diag_rhoi  ← from model
#   col 3: Nitot   col 7: diag_di    ← from model
_CAT_BLOCK_START = 15
_CAT_BLOCK_WIDTH =  8

# offsets within a per-category block
_OFF_QITOT = 0
_OFF_QIRIM = 1
_OFF_QILIQ = 2
_OFF_NITOT = 3
_OFF_BIRIM = 4
_OFF_ZITOT = 5
_OFF_RHOI  = 6   # model diag_rhoi
_OFF_DI    = 7   # model diag_di


def _cat_col(cat_idx: int, offset: int) -> int:
    """0-based column index for category cat_idx (0-based) and intra-block offset."""
    return _CAT_BLOCK_START + cat_idx * _CAT_BLOCK_WIDTH + offset


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def read_out_p3dat(filepath: str, n_cat: int = 1) -> np.ndarray:
    """
    Load out_p3.dat into a 2D numpy array of shape (n_levels, n_cols).

    Raises ValueError if the file has fewer columns than expected for n_cat.
    """
    data = np.loadtxt(filepath)
    if data.ndim == 1:
        data = data[np.newaxis, :]   # single-level edge case
    expected = _CAT_BLOCK_START + n_cat * _CAT_BLOCK_WIDTH
    if data.shape[1] < expected:
        raise ValueError(
            f'{filepath}: expected >= {expected} columns for nCat={n_cat}, '
            f'found {data.shape[1]}.')
    return data


def write_comparison_dat(
        filepath:   str,
        data:       np.ndarray,
        pp_results: dict
) -> None:
    """
    Write comparison output file in the same e13.4 format as out_p3.dat,
    with 6×n_cat post-processor columns appended after the original data.

    pp_results: dict returned by compute(), must contain keys
        'pp_rhoi', 'pp_di', 'pp_effi', 'pp_vmi', 'pp_lami', 'pp_mui'
        each with shape (n_levels, n_cat).

    Column layout of out_pp.dat (e.g. nCat=2):
        cols 0 – N-1 : all original columns from out_p3.dat  (N = data.shape[1])
        then per category (block of 6, cat 1 first):
          pp_rhoi  [kg m-3]
          pp_di    [m]
          pp_effi  [m]
          pp_vmi   [m s-1]   (rhofaci=1; valid at surface)
          pp_lami  [m-1]
          pp_mui   [-]
    """
    n_orig = data.shape[1]
    pp_rhoi = pp_results['pp_rhoi']
    pp_di   = pp_results['pp_di']
    pp_effi = pp_results['pp_effi']
    pp_vmi  = pp_results['pp_vmi']
    pp_lami = pp_results['pp_lami']
    pp_mui  = pp_results['pp_mui']
    n_cat   = pp_rhoi.shape[1]

    # interleave pp columns: 6 fields per category
    pp_cols = np.empty((data.shape[0], 6 * n_cat))
    for iice in range(n_cat):
        base = 6 * iice
        pp_cols[:, base + 0] = pp_rhoi[:, iice]
        pp_cols[:, base + 1] = pp_di[:,   iice]
        pp_cols[:, base + 2] = pp_effi[:, iice]
        pp_cols[:, base + 3] = pp_vmi[:,  iice]
        pp_cols[:, base + 4] = pp_lami[:, iice]
        pp_cols[:, base + 5] = pp_mui[:,  iice]

    out = np.column_stack([data, pp_cols])
    with open(filepath, 'w') as fh:
        for row in out:
            fh.write(''.join(f'{v:13.4e}' for v in row) + '\n')
    print(f'  Wrote {out.shape[0]} levels → {filepath}')
    for iice in range(n_cat):
        c = n_orig + 6 * iice
        print(f'  col {c}    = pp_rhoi  cat {iice+1}  [kg m-3]')
        print(f'  col {c+1}  = pp_di    cat {iice+1}  [m]')
        print(f'  col {c+2}  = pp_effi  cat {iice+1}  [m]')
        print(f'  col {c+3}  = pp_vmi   cat {iice+1}  [m s-1]  (rhofaci=1)')
        print(f'  col {c+4}  = pp_lami  cat {iice+1}  [m-1]')
        print(f'  col {c+5}  = pp_mui   cat {iice+1}  [-]')


# ---------------------------------------------------------------------------
# Python API
# ---------------------------------------------------------------------------
def compute(
        arrays:  dict,
        cfg:     P3Config,
        lut_dir: str  = 'lookup_tables',
        lut:     'P3LUT | None' = None
) -> dict:
    """
    Compute P3 ice diagnostics from prognostic arrays.

    This is the primary Python API entry point.  It accepts NumPy arrays
    directly, with no file I/O, and returns results as a dict of arrays.

    Parameters
    ----------
    arrays : dict
        Must contain the following keys (values are NumPy arrays):

        'Qitot'  [kg kg-1]  total ice mass mixing ratio
        'Qirim'  [kg kg-1]  rime ice mass mixing ratio
        'Qiliq'  [kg kg-1]  liquid-on-ice mass mixing ratio
        'Nitot'  [kg-1]     ice number concentration
        'Birim'  [m3 kg-1]  rime volume mixing ratio
        'Zitot'  [m6 kg-1]  reflectivity moment (required even if triple_moment=False)

        Each array may be:
          - 1-D, shape (n_levels,)          — any n_cat; interpreted as the
                                              *sum* / single-category field.
          - 2-D, shape (n_levels, n_cat)    — one column per ice category.

        For n_cat=1 a 1-D array is automatically reshaped to (n_levels, 1).
        For n_cat>1 a 2-D array of the correct shape is expected.

    cfg : P3Config
        Configuration (n_cat, triple_moment, liq_frac).

    lut_dir : str, optional
        Path to the directory containing unzipped LUT files.
        Ignored when *lut* is supplied.  Default: 'lookup_tables'.

    lut : P3LUT or None, optional
        Pre-loaded LUT instance.  Pass this to skip the (expensive) LUT
        loading step when calling compute() repeatedly.  When None (default)
        a new P3LUT is created from *lut_dir*.

    Returns
    -------
    dict with keys:
        'pp_rhoi'  ndarray, shape (n_levels, n_cat)  [kg m-3]
        'pp_di'    ndarray, shape (n_levels, n_cat)  [m]
        'pp_effi'  ndarray, shape (n_levels, n_cat)  [m]
        'pp_vmi'   ndarray, shape (n_levels, n_cat)  [m s-1]  (rhofaci=1; valid at surface)
        'pp_lami'  ndarray, shape (n_levels, n_cat)  [m-1]
        'pp_mui'   ndarray, shape (n_levels, n_cat)  [-]

    Raises
    ------
    ValueError
        If a required key is missing from *arrays*, or if an array has an
        unexpected shape.

    Examples
    --------
    >>> import numpy as np
    >>> from p3_postproc import compute, P3Config
    >>> cfg = P3Config(n_cat=1, triple_moment=False, liq_frac=False)
    >>> arrays = {k: np.zeros(50) for k in ('Qitot','Qirim','Qiliq','Nitot','Birim','Zitot')}
    >>> result = compute(arrays, cfg, lut_dir='lookup_tables')
    >>> result['pp_rhoi'].shape
    (50, 1)
    """
    # ---- input validation ----
    required = ('Qitot', 'Qirim', 'Qiliq', 'Nitot', 'Birim', 'Zitot')
    for key in required:
        if key not in arrays:
            raise ValueError(f"compute(): required key '{key}' missing from arrays dict.")

    # Normalise to 2-D (n_levels, n_cat)
    def _to_2d(arr, name):
        arr = np.asarray(arr, dtype=float)
        if arr.ndim == 1:
            arr = arr[:, np.newaxis]
        if arr.ndim != 2:
            raise ValueError(
                f"compute(): arrays['{name}'] must be 1-D or 2-D, got shape {arr.shape}.")
        return arr

    Qitot = _to_2d(arrays['Qitot'], 'Qitot')
    Qirim = _to_2d(arrays['Qirim'], 'Qirim')
    Qiliq = _to_2d(arrays['Qiliq'], 'Qiliq')
    Nitot = _to_2d(arrays['Nitot'], 'Nitot')
    Birim = _to_2d(arrays['Birim'], 'Birim')
    Zitot = _to_2d(arrays['Zitot'], 'Zitot')

    n_levels, n_cols = Qitot.shape
    if n_cols != cfg.n_cat:
        raise ValueError(
            f"compute(): arrays have {n_cols} column(s) but cfg.n_cat={cfg.n_cat}.")

    for name, arr in (('Qirim', Qirim), ('Qiliq', Qiliq),
                      ('Nitot', Nitot), ('Birim', Birim), ('Zitot', Zitot)):
        if arr.shape != (n_levels, cfg.n_cat):
            raise ValueError(
                f"compute(): arrays['{name}'].shape={arr.shape} does not match "
                f"Qitot.shape={Qitot.shape}.")

    # If liqFrac is disabled, treat Qiliq as zero (avoids caller needing to zero it)
    if not cfg.liq_frac:
        Qiliq = np.zeros_like(Qiliq)

    # ---- LUT initialisation ----
    if lut is None:
        lut = P3LUT(lut_dir, cfg)

    # ---- compute diagnostics ----
    pp_rhoi = np.zeros((n_levels, cfg.n_cat))
    pp_di   = np.zeros((n_levels, cfg.n_cat))
    pp_effi = np.zeros((n_levels, cfg.n_cat))
    pp_vmi  = np.zeros((n_levels, cfg.n_cat))
    pp_lami = np.zeros((n_levels, cfg.n_cat))
    pp_mui  = np.zeros((n_levels, cfg.n_cat))

    for k in range(n_levels):
        for iice in range(cfg.n_cat):
            (pp_rhoi[k, iice], pp_di[k, iice],
             pp_effi[k, iice], pp_vmi[k, iice],
             pp_lami[k, iice], pp_mui[k, iice]) = lut.get_diag(
                Qitot[k, iice], Qirim[k, iice], Qiliq[k, iice],
                Nitot[k, iice], Birim[k, iice], Zitot[k, iice])

    return {
        'pp_rhoi': pp_rhoi,
        'pp_di':   pp_di,
        'pp_effi': pp_effi,
        'pp_vmi':  pp_vmi,
        'pp_lami': pp_lami,
        'pp_mui':  pp_mui,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='P3 offline post-processor — compute diag_rhoi and diag_di '
                    'from prognostic ice variables and compare to model output.')
    parser.add_argument('--input',   default='out_p3.dat',
                        help='Input file from cld1d   (default: out_p3.dat)')
    parser.add_argument('--lut-dir', default='lookup_tables',
                        help='Directory with unzipped LUT files  (default: lookup_tables)')
    parser.add_argument('--config',  default='p3_config.txt',
                        help='P3 config file           (default: p3_config.txt)')
    parser.add_argument('--output',  default='out_pp.dat',
                        help='Output comparison file   (default: out_pp.dat)')
    args = parser.parse_args()

    print('\nP3 Post-Processor  (p3_postproc.py)')
    print(f'  input   : {args.input}')
    print(f'  lut-dir : {args.lut_dir}')
    print(f'  config  : {args.config}')
    print(f'  output  : {args.output}')

    # --- configuration ---
    cfg = P3Config.from_file(args.config)
    print(f'  nCat={cfg.n_cat}  3-moment={cfg.triple_moment}  liqFrac={cfg.liq_frac}')

    # --- load LUT ---
    lut = P3LUT(args.lut_dir, cfg)

    # --- read input data ---
    print(f'  Reading {args.input} ...')
    data     = read_out_p3dat(args.input, cfg.n_cat)
    n_levels = data.shape[0]
    print(f'  {n_levels} levels, {data.shape[1]} columns')

    # --- build per-category arrays and call compute() ---
    arrays = {}
    for name, off in (('Qitot', _OFF_QITOT), ('Qirim', _OFF_QIRIM),
                      ('Qiliq', _OFF_QILIQ), ('Nitot', _OFF_NITOT),
                      ('Birim', _OFF_BIRIM), ('Zitot', _OFF_ZITOT)):
        arrays[name] = np.column_stack(
            [data[:, _cat_col(iice, off)] for iice in range(cfg.n_cat)])

    result  = compute(arrays, cfg, lut=lut)
    pp_rhoi = result['pp_rhoi']
    pp_di   = result['pp_di']

    n_active = np.zeros(cfg.n_cat, dtype=int)
    for iice in range(cfg.n_cat):
        n_active[iice] = int(
            np.sum(data[:, _cat_col(iice, _OFF_QITOT)] >= QSMALL))
        print(f'  Post-processor: {n_active[iice]}/{n_levels} active levels '
              f'(cat {iice+1})')

    # --- write comparison output ---
    write_comparison_dat(args.output, data, result)

    # --- verification statistics (per category) ---
    print(f'\n  Verification statistics:')
    print(f'    {"":16s}  {"max |diff|":>12s}  {"mean |diff|":>12s}  '
          f'{"max |rel err|":>14s}  {"median |rel err|":>16s}')

    # Column indices for extended model reference fields (nCat=1 only, when present).
    # Written by cld1d when the build includes the verification patch.
    _n_std_cols = _CAT_BLOCK_START + cfg.n_cat * _CAT_BLOCK_WIDTH  # 23 for nCat=1
    _has_ext = (cfg.n_cat == 1 and data.shape[1] >= _n_std_cols + 4)
    _COL_MODEL_EFFI = _n_std_cols + 0   # col 23
    _COL_MODEL_VMI  = _n_std_cols + 1   # col 24
    _COL_MODEL_LAMI = _n_std_cols + 2   # col 25
    _COL_MODEL_MUI  = _n_std_cols + 3   # col 26

    def _stats(diff, ref):
        abs_rel = np.abs(diff / np.where(ref != 0, ref, np.nan))
        mx      = np.max(np.abs(diff))
        mn      = np.mean(np.abs(diff))
        rel_max = np.nanmax(abs_rel)
        rel_med = np.nanmedian(abs_rel)
        return mx, mn, rel_max, rel_med

    any_active = False
    for iice in range(cfg.n_cat):
        active = data[:, _cat_col(iice, _OFF_QITOT)] >= QSMALL
        if not active.any():
            continue
        any_active = True
        cat_label  = f'cat {iice+1}' if cfg.n_cat > 1 else ''

        model_rhoi = data[active, _cat_col(iice, _OFF_RHOI)]
        model_di   = data[active, _cat_col(iice, _OFF_DI)]
        mx, mn, rel_max, rel_med = _stats(pp_rhoi[active, iice] - model_rhoi, model_rhoi)
        print(f'    {f"diag_rhoi {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel_max:14.4e}  {rel_med:16.4e}')
        mx, mn, rel_max, rel_med = _stats(pp_di[active, iice] - model_di, model_di)
        print(f'    {f"diag_di   {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel_max:14.4e}  {rel_med:16.4e}')

        if _has_ext and iice == 0:
            model_effi = data[active, _COL_MODEL_EFFI]
            model_vmi  = data[active, _COL_MODEL_VMI]
            model_lami = data[active, _COL_MODEL_LAMI]
            model_mui  = data[active, _COL_MODEL_MUI]
            mx, mn, rel_max, rel_med = _stats(result['pp_effi'][active, iice] - model_effi, model_effi)
            print(f'    {f"diag_effi {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel_max:14.4e}  {rel_med:16.4e}')
            mx, mn, rel_max, rel_med = _stats(result['pp_vmi'][active, iice]  - model_vmi,  model_vmi)
            print(f'    {f"diag_vmi  {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel_max:14.4e}  {rel_med:16.4e}  (rhofaci=1 vs model)')
            mx, mn, rel_max, rel_med = _stats(result['pp_lami'][active, iice] - model_lami, model_lami)
            print(f'    {f"diag_lami {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel_max:14.4e}  {rel_med:16.4e}')
            mx, mn, rel_max, rel_med = _stats(result['pp_mui'][active, iice]  - model_mui,  model_mui)
            print(f'    {f"diag_mui  {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel_max:14.4e}  {rel_med:16.4e}')

    if not any_active:
        print('  No active ice levels found — nothing to compare.')


if __name__ == '__main__':
    main()
