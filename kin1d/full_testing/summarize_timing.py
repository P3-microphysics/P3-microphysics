#!/usr/bin/env python3
"""
summarize_timing.py — parse and summarize kin1d multi-run timing results.

Usage
-----
  python summarize_timing.py [--dir DIR] [--hash HASH]

Options
-------
  --dir  DIR   Directory containing timing_<config>_run<N>.txt files
               [default: kin1d_tests/timing-results]
  --hash HASH  Git hash label for the report header [default: auto-detect]
  --save PATH  Write report to PATH in addition to stdout
"""

import argparse
import glob
import os
import re
import sys

CONFIGS = ['1FF', '1FT', '1TF', '1TT', '2TT', '3TT']

BLOCK_LABELS = {
    1: 'full p3_main',
    2: 'i_loop_main',
    3: 'k_loop_main_process',
    6: 'sedimentation',
    9: 'type_diags',
}


def parse_timing_file(path):
    """Return dict of {block_num: time_ms} from a timing.txt file."""
    result = {}
    in_block = False
    with open(path) as f:
        for line in f:
            if 'Block Description' in line:
                in_block = True
                continue
            if in_block and line.strip().startswith('---'):
                continue
            if in_block and line.strip() == '':
                break
            if in_block:
                # e.g. "   1 full p3_main              11.772900      100.000000"
                m = re.match(r'\s+(\d+)\s+.+?\s+([\d.]+)\s+([\d.]+)\s*$', line)
                if m:
                    block = int(m.group(1))
                    ms = float(m.group(2))
                    result[block] = ms
    return result


def collect(directory):
    """
    Return nested dict: data[config][block_num] = [ms_run1, ms_run2, ...]
    """
    data = {cfg: {} for cfg in CONFIGS}
    for cfg in CONFIGS:
        pattern = os.path.join(directory, f'timing_{cfg}_run*.txt')
        files = sorted(glob.glob(pattern))
        if not files:
            continue
        for path in files:
            parsed = parse_timing_file(path)
            for block, ms in parsed.items():
                data[cfg].setdefault(block, []).append(ms)
    return data


def stats(values):
    mean = sum(values) / len(values)
    return mean, min(values), max(values)


def format_report(data, git_hash, nruns):
    lines = []
    sep = '=' * 72
    lines.append(sep)
    lines.append(f'  P3 kin1d timing report')
    lines.append(f'  hash: {git_hash}   build: -O2   runs per config: {nruns}')
    lines.append(sep)

    # Headline table: Block 1 (full p3_main) across configs
    lines.append('')
    lines.append('  Full p3_main (Block 1) — ms')
    run_header = '  '.join([f'r{i+1:>4}' for i in range(nruns)])
    lines.append(f"  {'Config':<8}  {run_header}  {'Mean':>8}  {'Min':>8}  {'Max':>8}  {'Range':>8}")
    lines.append(f"  {'-'*8}  {'  '.join(['-'*5]*nruns)}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    for cfg in CONFIGS:
        vals = data[cfg].get(1)
        if not vals:
            lines.append(f"  {cfg:<8}  (no data)")
            continue
        mean, lo, hi = stats(vals)
        run_strs = '  '.join([f'{v:>5.2f}' for v in vals])
        lines.append(f"  {cfg:<8}  {run_strs}  {mean:>8.3f}  {lo:>8.3f}  {hi:>8.3f}  {hi-lo:>8.3f}")

    # Per-block detail for each config
    for cfg in CONFIGS:
        if not data[cfg]:
            continue
        lines.append('')
        lines.append(f'  Per-block detail — {cfg}  (mean across {nruns} runs, ms)')
        lines.append(f"  {'Block':<4}  {'Description':<26}  {'Mean':>8}  {'Min':>8}  {'Max':>8}  {'% of B1':>8}")
        lines.append(f"  {'-'*4}  {'-'*26}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")
        b1_mean = stats(data[cfg].get(1, [1]))[0]
        for block in sorted(data[cfg].keys()):
            vals = data[cfg][block]
            mean, lo, hi = stats(vals)
            label = BLOCK_LABELS.get(block, f'block {block}')
            pct = 100. * mean / b1_mean if b1_mean else 0.
            lines.append(f"  {block:<4}  {label:<26}  {mean:>8.3f}  {lo:>8.3f}  {hi:>8.3f}  {pct:>7.1f}%")

    lines.append('')
    lines.append(sep)
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Summarize kin1d multi-run timing results.')
    parser.add_argument('--dir',  default='kin1d_tests/timing-results')
    parser.add_argument('--hash', default=None)
    parser.add_argument('--save', default=None)
    args = parser.parse_args()

    if args.hash is None:
        try:
            import subprocess
            args.hash = subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'], stderr=subprocess.DEVNULL
            ).decode().strip()
        except Exception:
            args.hash = 'unknown'

    data = collect(args.dir)
    all_nruns = [len(v) for cfg in data.values() for v in cfg.values() if v]
    nruns = max(all_nruns) if all_nruns else 0

    if nruns == 0:
        print(f"No timing files found in {args.dir}", file=sys.stderr)
        sys.exit(1)

    report = format_report(data, args.hash, nruns)
    print(report)

    if args.save:
        with open(args.save, 'w') as f:
            f.write(report + '\n')
        print(f"\n  Report saved to: {args.save}")


if __name__ == '__main__':
    main()
