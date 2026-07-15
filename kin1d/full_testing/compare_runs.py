#!/usr/bin/env python3
"""
compare_runs.py — kin1d test output comparison tool for P3 microphysics.

Compares two sets of kin1d output files (out_p3_<ID>.dat) and reports
whether results are bit-identical, within tolerance, or differ significantly.

Usage
-----
  python compare_runs.py [options]

Options
-------
  --ref   DIR   Reference results directory  [default: kin1d_tests/results-reference]
  --new   DIR   New results directory        [default: kin1d_tests/results-last_commit]
  --mode  STR   Comparison mode:
                  exact  — require bit-identical match (default)
                  tol    — allow differences up to --atol / --rtol
  --atol  FLOAT Absolute tolerance for 'tol' mode  [default: 1e-10]
  --rtol  FLOAT Relative tolerance for 'tol' mode  [default: 1e-5]
  --configs STR Comma-separated list of configs to compare [default: 1FF,1FT,1TF,1TT,2TT,3TT]

Exit code
---------
  0  All configs PASS
  1  One or more configs FAIL

Column layout (from cld1d.f90)
-------------------------------
  0  z        height AGL (m)
  1  w        vertical velocity (m/s)
  2  prt_liq  liquid precip rate (mm/h)
  3  prt_sol  solid precip rate (mm/h)
  4  Ze       equivalent reflectivity (dBZ)
  5  T-T0     temperature relative to 0C (K)
  6  Qc       cloud water mixing ratio
  7  Qr       rain mixing ratio
  8  Nc       cloud droplet number
  9  Nr       rain number
 10  Qi_tot   total ice mixing ratio (all cats summed)
 11  Ni_tot   total ice number (all cats summed)
 12  Frime    bulk rime fraction (cat 1)
 13  Fliq     bulk liquid fraction (cat 1)
 14  Drm      mean-mass rain diameter
 15  Qitot_1  ice mass (cat 1)
 16  Qirim_1  rime mass (cat 1)
 17  Qiliq_1  liquid mass on ice (cat 1)
 18  Nitot_1  ice number (cat 1)
 19  Birim_1  rime volume (cat 1)
 20  Zitot_1  6th moment (cat 1)
 21  rhoi_1   bulk ice density (cat 1)
 22  di_1     mean ice diameter (cat 1)
 [nCat>=2: cols 23-30 repeat for cat 2]
 [nCat>=3: cols 31-38 repeat for cat 3]
"""

import argparse
import os
import sys
import numpy as np

# ── column names (shared prefix for per-category ice variables) ──────────────
COL_NAMES_BASE = [
    'z', 'w', 'prt_liq', 'prt_sol', 'Ze', 'T-T0',
    'Qc', 'Qr', 'Nc', 'Nr', 'Qi_sum', 'Ni_sum',
    'Frime', 'Fliq', 'Drm',
]
COL_NAMES_ICE = ['Qitot', 'Qirim', 'Qiliq', 'Nitot', 'Birim', 'Zitot', 'rhoi', 'di']

# Fill-value used for "sentinel" Ze = -99 (no precipitation)
ZE_SENTINEL = -99.0


def col_names_for_ncols(ncols):
    """Return a list of column name strings for a given number of output columns."""
    names = list(COL_NAMES_BASE)
    col = len(names)
    cat = 1
    while col < ncols:
        for var in COL_NAMES_ICE:
            names.append(f'{var}_{cat}')
            col += 1
            if col >= ncols:
                break
        cat += 1
    return names


def load(path):
    """Load an out_p3.dat file; return a 2-D numpy array (rows × cols)."""
    try:
        data = np.loadtxt(path)
    except Exception as e:
        print(f"  ERROR: cannot read {path}: {e}")
        return None
    if data.ndim == 1:
        print(f"  ERROR: {path} appears empty or has only one row.")
        return None
    return data


def compare_exact(ref, new, col_names, config_id):
    """Exact (bit-identical) comparison. Returns (passed, details_str)."""
    if ref.shape != new.shape:
        return False, f"  Shape mismatch: ref={ref.shape}  new={new.shape}"

    match = np.array_equal(ref, new)
    if match:
        return True, "  All values bit-identical ✓"

    # Find where they differ
    diff_mask = ref != new
    n_diff = np.sum(diff_mask)
    rows_diff, cols_diff = np.where(diff_mask)
    unique_cols = sorted(set(cols_diff))
    col_labels = [col_names[c] if c < len(col_names) else f'col{c}' for c in unique_cols]

    abs_diff = np.abs(ref - new)
    max_abs = np.max(abs_diff)
    # relative diff — avoid div-by-zero with small reference values
    with np.errstate(invalid='ignore', divide='ignore'):
        rel_diff = np.where(np.abs(ref) > 1e-30,
                            abs_diff / np.abs(ref),
                            abs_diff)
    max_rel = np.max(rel_diff)

    lines = [
        f"  FAIL — {n_diff} of {ref.size} values differ",
        f"  max |Δ| = {max_abs:.3e}   max |Δ/ref| = {max_rel:.3e}",
        f"  Columns with differences: {', '.join(col_labels)}",
    ]

    # Per-column breakdown for differing columns
    for c in unique_cols[:8]:   # cap at 8 to keep output readable
        label = col_names[c] if c < len(col_names) else f'col{c}'
        col_abs = np.max(abs_diff[:, c])
        with np.errstate(invalid='ignore', divide='ignore'):
            col_rel = np.max(rel_diff[:, c])
        lines.append(f"    {label:15s}  max|Δ|={col_abs:.3e}  max|Δ/ref|={col_rel:.3e}")
    if len(unique_cols) > 8:
        lines.append(f"    ... and {len(unique_cols)-8} more columns")

    return False, '\n'.join(lines)


def compare_tol(ref, new, col_names, config_id, atol, rtol):
    """Tolerance-based comparison. Returns (passed, details_str)."""
    if ref.shape != new.shape:
        return False, f"  Shape mismatch: ref={ref.shape}  new={new.shape}"

    abs_diff = np.abs(ref - new)
    with np.errstate(invalid='ignore', divide='ignore'):
        rel_diff = np.where(np.abs(ref) > 1e-30,
                            abs_diff / np.abs(ref),
                            abs_diff)

    fail_mask = (abs_diff > atol) & (rel_diff > rtol)
    n_fail = np.sum(fail_mask)

    max_abs_global = float(np.max(abs_diff))
    max_rel_global = float(np.max(rel_diff))

    lines = [
        f"  max |Δ|     = {max_abs_global:.3e}  (atol={atol:.1e})",
        f"  max |Δ/ref| = {max_rel_global:.3e}  (rtol={rtol:.1e})",
    ]

    # Per-column summary
    ncols = ref.shape[1]
    col_stats = []
    for c in range(ncols):
        label = col_names[c] if c < len(col_names) else f'col{c}'
        col_abs = float(np.max(abs_diff[:, c]))
        col_rel = float(np.max(rel_diff[:, c]))
        col_fail = int(np.sum(fail_mask[:, c]))
        if col_abs > 0:
            col_stats.append((col_abs, col_rel, col_fail, label, c))

    col_stats.sort(reverse=True)
    if col_stats:
        lines.append(f"  {'Column':<16}  {'max|Δ|':>12}  {'max|Δ/ref|':>12}  {'#outside_tol':>13}")
        lines.append(f"  {'-'*16}  {'-'*12}  {'-'*12}  {'-'*13}")
        for col_abs, col_rel, col_fail, label, _ in col_stats[:12]:
            flag = '  *** ' if col_fail > 0 else ''
            lines.append(f"  {label:<16}  {col_abs:>12.3e}  {col_rel:>12.3e}  {col_fail:>13d}{flag}")
        if len(col_stats) > 12:
            lines.append(f"  ... ({len(col_stats)-12} more columns with non-zero diffs)")

    passed = n_fail == 0
    status = "PASS" if passed else f"FAIL ({n_fail} values outside tolerance)"
    lines.insert(0, f"  {status}")
    return passed, '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Compare two sets of kin1d output files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--ref',     default='kin1d_tests/results-reference',
                        help='Reference results directory')
    parser.add_argument('--new',     default='kin1d_tests/results-last_commit',
                        help='New results directory')
    parser.add_argument('--mode',    default='exact', choices=['exact', 'tol'],
                        help='Comparison mode: exact (bit-match) or tol (tolerance)')
    parser.add_argument('--atol',    type=float, default=1e-10,
                        help='Absolute tolerance [tol mode only]')
    parser.add_argument('--rtol',    type=float, default=1e-5,
                        help='Relative tolerance [tol mode only]')
    parser.add_argument('--configs', default='1FF,1FT,1TF,1TT,2TT,3TT',
                        help='Comma-separated config IDs to compare')
    args = parser.parse_args()

    configs = [c.strip() for c in args.configs.split(',')]

    print()
    print('=' * 68)
    print(f'  P3 kin1d results comparison   mode={args.mode}')
    print(f'  ref: {args.ref}')
    print(f'  new: {args.new}')
    if args.mode == 'tol':
        print(f'  atol={args.atol:.1e}   rtol={args.rtol:.1e}')
    print('=' * 68)

    results = {}
    for cfg in configs:
        ref_path = os.path.join(args.ref, f'out_p3_{cfg}.dat')
        new_path = os.path.join(args.new, f'out_p3_{cfg}.dat')

        print(f'\n  Config {cfg}')
        print(f'  {"-" * 40}')

        if not os.path.exists(ref_path):
            print(f'  SKIP — reference file not found: {ref_path}')
            results[cfg] = None
            continue
        if not os.path.exists(new_path):
            print(f'  SKIP — new file not found: {new_path}')
            results[cfg] = None
            continue

        ref_data = load(ref_path)
        new_data = load(new_path)
        if ref_data is None or new_data is None:
            results[cfg] = False
            continue

        col_names = col_names_for_ncols(ref_data.shape[1])

        if args.mode == 'exact':
            passed, details = compare_exact(ref_data, new_data, col_names, cfg)
        else:
            passed, details = compare_tol(ref_data, new_data, col_names, cfg,
                                          atol=args.atol, rtol=args.rtol)

        print(details)
        results[cfg] = passed

    # Summary
    print()
    print('=' * 68)
    print('  SUMMARY')
    print('=' * 68)
    all_passed = True
    for cfg in configs:
        r = results.get(cfg)
        if r is True:
            status = 'PASS'
        elif r is False:
            status = 'FAIL'
            all_passed = False
        else:
            status = 'SKIP'
        print(f'  {cfg:8s}  {status}')

    overall = 'ALL PASS' if all_passed else 'FAILURES DETECTED'
    print()
    print(f'  Result: {overall}')
    print('=' * 68)
    print()

    sys.exit(0 if all_passed else 1)


if __name__ == '__main__':
    main()
