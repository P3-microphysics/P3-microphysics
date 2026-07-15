#!/bin/bash
# save_reference.sh — save the current kin1d test output as the reference baseline.
#
# Run this BEFORE making code changes to establish the baseline you want to
# compare against.  The reference is stored in kin1d_tests/results-reference/.
#
# Usage (from kin1d/src/):
#   ./save_reference.sh              # run tests and save as reference
#   ./save_reference.sh --no-run     # just copy last_commit -> reference
#                                    #   (if you already ran run_full_tests)
#
# After making your code changes, run:
#   ./run_full_tests --compare       # for changes expected to be bit-identical
#   ./run_full_tests --compare-tol   # for changes expected to change results slightly

set -e

rootdir="."
NO_RUN=0
for arg in "$@"; do
    [ "$arg" = "--no-run" ] && NO_RUN=1
done

if [ "$NO_RUN" = "0" ]; then
    echo "Running full test suite to generate reference outputs..."
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    bash "${SCRIPT_DIR}/run_full_tests"
fi

SRC="${rootdir}/kin1d_tests/results-last_commit"
DST="${rootdir}/kin1d_tests/results-reference"

if [ ! -d "$SRC" ]; then
    echo "ERROR: $SRC not found. Run run_full_tests first."
    exit 1
fi

rm -rf "$DST"
cp -r "$SRC" "$DST"

echo
echo "Reference baseline saved to: $DST"
echo "Contents:"
ls "$DST"
echo
echo "This baseline will be used by compare_runs.py."
echo "Run 'run_full_tests --compare' after making changes."
