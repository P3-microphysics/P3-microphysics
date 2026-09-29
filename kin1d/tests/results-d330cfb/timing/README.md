# Timing: 5.6.0+qwgrth (d330cfb) vs v5.6.0 (dev, src identical to 390c9d6)

kin1d, `-O2 -cpp -DTIMING_P3`, interleaved A/B (ref, cand, ref, cand, ...), 8 runs per build per
config, Mac (Apple Silicon, x86_64 gfortran 14.2 under Rosetta), machine otherwise idle
(load ~2.5). Both builds used the same kin1d scripts, configs and Makefile (from bugfix-qwgrth).
Values: block 1 (full `p3_main`) as reported in `timing.txt`.

| Config | ref mean (SD) | cand mean (SD) | Δ mean | Δ min |
|---|---|---|---|---|
| 1FF | 0.4379 (5.4 %) | 0.4292 (5.2 %) | −2.0 % | −0.4 % |
| 1FT | 0.4504 (4.2 %) | 0.4546 (4.5 %) | +0.9 % | +1.8 % |
| 1TF | 0.4952 (2.6 %) | 0.4912 (7.2 %) | −0.8 % | −0.8 % |
| 1TT | 0.5021 (6.3 %) | 0.5037 (4.5 %) | +0.3 % | +0.4 % |
| 2TT | 0.5842 (2.4 %) | 0.6171 (2.1 %) | +5.6 % | +3.9 % |
| 3TT | 0.8191 (3.1 %) | 0.8456 (2.5 %) | +3.2 % | +5.1 % |

Interpretation: single-category configs are timing-neutral (within run-to-run noise). The
multi-category configs are 3–6 % slower, consistently in mean and minimum. The code change adds
only a few operations per wet-growth point; the multi-category kin1d solutions differ strongly
with the fix (2TT 90-min precipitation 3.2 → 9.1 mm), so more points carry active rain/ice and
the extra cost most likely reflects the changed physical state rather than the formulas. This
test cannot separate the two.

Files: `ref-dev/`, `cand/` (per-run timing files), `timing_report_*.txt` (summarize_timing.py).
