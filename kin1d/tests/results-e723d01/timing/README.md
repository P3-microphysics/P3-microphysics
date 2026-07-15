# Timing comparison: 4-routine sedimentation (77463c5) vs merged sedimentation_ice (e723d01)

Both reports come from an interleaved campaign (-O2, no bounds-check): for each of
the 6 configs, 8 alternating runs of the reference (77463c5, four sedimentation_ice_*
routines) and merged (e723d01, single flag-driven sedimentation_ice) executables,
so both builds sample identical machine conditions.

Summary (total p3_main, mean of 8 runs, ms):

| Config | 77463c5 | e723d01 | delta  |
|--------|---------|---------|--------|
| 1FF    | 0.421   | 0.434   | +3.1%* |
| 1FT    | 0.442   | 0.453   | +2.5%  |
| 1TF    | 0.493   | 0.492   | -0.2%  |
| 1TT    | 0.501   | 0.499   | -0.4%  |
| 2TT    | 0.610   | 0.602   | -1.3%  |
| 3TT    | 0.825   | 0.829   | +0.5%  |

*one 0.495 outlier; min-to-min +2.4%.

All deltas are smaller than one within-build run-to-run SD (2-5%). The instrumented
sedimentation block is a consistent +~2 us/timestep (~3% of the block). Per-config
minima bound the true residual cost at <=~2% of p3_main.

NOTE: timing still to be verified in 3D models (GEM/WRF). If a regression appears
there, Option B (include-file compile-time specialization) is the documented
fallback — see CLAUDE.md "Sedimentation Routine Consolidation".
