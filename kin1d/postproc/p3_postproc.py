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

The output file has the same format as out_p3.dat with two extra columns:
    col N  : pp_rhoi  [kg m-3]  ← post-processor bulk ice density
    col N+1: pp_di    [m]       ← post-processor mean-mass diameter

The original model diag_rhoi and diag_di columns are preserved, allowing
direct side-by-side comparison. Expected accuracy: median relative error ~0.03%.

COMMAND-LINE OPTIONS:

    --input FILE           Input file (default: out_p3.dat)
    --lut-dir DIR          LUT directory (default: lookup_tables)
    --config FILE          P3 config file (default: p3_config.txt)
    --output FILE          Output file (default: out_pp.dat)

SUPPORTED CONFIGURATIONS (v0.2):

    nCat=1–3, trplMomIce=.false., liqFrac=.false.  (2-moment, no liquid)
    nCat=1–3, trplMomIce=.false., liqFrac=.true.   (2-moment, with liquid fraction)

3-moment ice (trplMomIce=.true.) is not yet implemented.

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

# 2-moment LUT dimensions  (integer, parameter :: ... in microphy_p3.f90)
ISIZE     = 50    # particle size bins
RIMSIZE   =  4    # rime mass fraction bins
LIQSIZE   =  4    # liquid mass fraction bins
DENSIZE   =  5    # bulk rime density bins
TABSIZE   = 19    # number of process quantities in 2-moment table
RCOLLSIZE = 30    # ice-rain collection bins (file structure only; not used here)

# Process indices in the LUT (1-based, matching Fortran)
IDX_DIAG_DI   = 11    # mean-mass diameter    [m]
IDX_DIAG_RHOI = 12    # bulk ice density      [kg m-3]


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

    Implements the 2-moment ice path, mirroring:
        calc_bulkRhoRime()          in microphy_p3.f90
        find_lookupTable_indices_1a()  "
        proc_from_LUT_main2mom()       "
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
            lut_file = candidates[0]
            self._load_2mom_lut(lut_file)
        else:
            raise NotImplementedError('3-moment ice LUT not yet implemented.')

    # ------------------------------------------------------------------
    def _load_2mom_lut(self, path: Path):
        """
        Read the 2-moment LUT into self.itab[densize, rimsize, liqsize, isize, tabsize].

        File structure mirrors the Fortran p3_init read loop:
          - 1 header line (skipped)
          - For each (jj, ii, ll):
              isize rows × 25 cols : 4 dummy | 19 table values | 2 dummy
              isize×rcollsize rows  : ice-rain collection data (skipped)
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
    def get_diag(
            self,
            qi_tot: float, qi_rim: float, qi_liq: float,
            ni_tot: float, bi_rim: float
    ) -> tuple[float, float]:
        """
        Compute diag_rhoi and diag_di for a single grid point.

        Replicates the qi_not_small block in p3_main (microphy_p3.f90).

        Returns (0., 0.) when qi_tot < QSMALL, matching model behaviour.
        """
        if qi_tot < QSMALL:
            return 0., 0.

        ni_tot = max(ni_tot, NSMALL)

        # calc_bulkRhoRime may adjust qi_rim and bi_rim; use returned values
        rho_p, qi_rim, bi_rim = self.calc_bulk_rho_rime(
            qi_tot, qi_rim, qi_liq, bi_rim)

        dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7 = self._find_indices(
            qi_tot, ni_tot, qi_rim, qi_liq, rho_p)

        diag_rhoi = self._interp_2mom(
            IDX_DIAG_RHOI, dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7)
        diag_di   = self._interp_2mom(
            IDX_DIAG_DI,   dumi, dumjj, dumii, dumll, dum1, dum5, dum4, dum7)

        return diag_rhoi, diag_di


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
        filepath: str,
        data: np.ndarray,
        pp_rhoi: np.ndarray,
        pp_di:   np.ndarray
) -> None:
    """
    Write comparison output file in the same e13.4 format as out_p3.dat,
    with 2×n_cat post-processor columns appended after the original data.

    pp_rhoi, pp_di: shape (n_levels, n_cat)

    Column layout of out_pp.dat (e.g. nCat=2):
        cols 0 – N-1 : all original columns from out_p3.dat  (N = data.shape[1])
        col  N       : pp_rhoi  cat 1  [kg m-3]
        col  N+1     : pp_di    cat 1  [m]
        col  N+2     : pp_rhoi  cat 2  [kg m-3]
        col  N+3     : pp_di    cat 2  [m]
        ...
    """
    n_orig = data.shape[1]
    n_cat  = pp_rhoi.shape[1]

    # interleave pp columns: rhoi_1, di_1, rhoi_2, di_2, ...
    pp_cols = np.empty((data.shape[0], 2 * n_cat))
    for iice in range(n_cat):
        pp_cols[:, 2 * iice]     = pp_rhoi[:, iice]
        pp_cols[:, 2 * iice + 1] = pp_di[:,   iice]

    out = np.column_stack([data, pp_cols])
    with open(filepath, 'w') as fh:
        for row in out:
            fh.write(''.join(f'{v:13.4e}' for v in row) + '\n')
    print(f'  Wrote {out.shape[0]} levels → {filepath}')
    for iice in range(n_cat):
        c = n_orig + 2 * iice
        print(f'  col {c}  = pp_rhoi  cat {iice+1}  [kg m-3]')
        print(f'  col {c+1}  = pp_di    cat {iice+1}  [m]')


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

    if cfg.triple_moment:
        raise NotImplementedError('trplMomIce=.true.: 3-moment ice not yet implemented.')

    # --- load LUT ---
    lut = P3LUT(args.lut_dir, cfg)

    # --- read input data ---
    print(f'  Reading {args.input} ...')
    data     = read_out_p3dat(args.input, cfg.n_cat)
    n_levels = data.shape[0]
    print(f'  {n_levels} levels, {data.shape[1]} columns')

    # --- compute post-processor diagnostics at each level and category ---
    pp_rhoi  = np.zeros((n_levels, cfg.n_cat))
    pp_di    = np.zeros((n_levels, cfg.n_cat))
    n_active = np.zeros(cfg.n_cat, dtype=int)

    for k in range(n_levels):
        for iice in range(cfg.n_cat):
            qi_tot = data[k, _cat_col(iice, _OFF_QITOT)]
            qi_rim = data[k, _cat_col(iice, _OFF_QIRIM)]
            qi_liq = data[k, _cat_col(iice, _OFF_QILIQ)]
            ni_tot = data[k, _cat_col(iice, _OFF_NITOT)]
            bi_rim = data[k, _cat_col(iice, _OFF_BIRIM)]

            pp_rhoi[k, iice], pp_di[k, iice] = lut.get_diag(
                qi_tot, qi_rim, qi_liq, ni_tot, bi_rim)
            if qi_tot >= QSMALL:
                n_active[iice] += 1

    for iice in range(cfg.n_cat):
        print(f'  Post-processor: {n_active[iice]}/{n_levels} active levels '
              f'(cat {iice+1})')

    # --- write comparison output ---
    write_comparison_dat(args.output, data, pp_rhoi, pp_di)

    # --- verification statistics (per category) ---
    print(f'\n  Verification statistics:')
    print(f'    {"":16s}  {"max |diff|":>12s}  {"mean |diff|":>12s}  '
          f'{"max |rel err|":>14s}')

    any_active = False
    for iice in range(cfg.n_cat):
        active = data[:, _cat_col(iice, _OFF_QITOT)] >= QSMALL
        if not active.any():
            continue
        any_active = True
        model_rhoi = data[active, _cat_col(iice, _OFF_RHOI)]
        model_di   = data[active, _cat_col(iice, _OFF_DI)]
        diff_rhoi  = pp_rhoi[active, iice] - model_rhoi
        diff_di    = pp_di[active,   iice] - model_di
        cat_label  = f'cat {iice+1}' if cfg.n_cat > 1 else ''

        def _stats(diff, ref):
            mx  = np.max(np.abs(diff))
            mn  = np.mean(np.abs(diff))
            rel = np.max(np.abs(diff / np.where(ref != 0, ref, np.nan)))
            return mx, mn, rel

        mx, mn, rel = _stats(diff_rhoi, model_rhoi)
        print(f'    {f"diag_rhoi {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel:14.4e}')
        mx, mn, rel = _stats(diff_di, model_di)
        print(f'    {f"diag_di   {cat_label}":16s}  {mx:12.4e}  {mn:12.4e}  {rel:14.4e}')

    if not any_active:
        print('  No active ice levels found — nothing to compare.')


if __name__ == '__main__':
    main()
