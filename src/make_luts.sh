#!/bin/bash
# make_luts.sh -- generate the P3 lookup tables, with N blocks running in parallel.
#
#   make_luts.sh TABLE [-j N] [--fc FC] [--outdir DIR] [--keep]
#
#   TABLE     lt1-2mom | lt1-3mom | lt2 | lt3 | all
#   -j N      number of blocks run at once (default 1).  On a machine of your own, the number of
#             cores (-j $(getconf _NPROCESSORS_ONLN)); in a batch job, the cores requested; on a
#             shared login node, never more than 1-2.
#   --fc FC   ifx | ifort | gfortran (default: the first one found, in that order)
#   --outdir  where the table(s) and their .info files go (default: current directory)
#   --keep    keep the work directory (blocks, logs) after a successful run
#   -h        this help
#
# Examples:
#   src/make_luts.sh lt3 -j 8 --outdir ~/p3-lut/test            # minutes on 8 cores
#   src/make_luts.sh all -j 128 --fc ifx --outdir /some/dir     # inside a 128-core batch job
#
# What it does.  Each generator (src/create_p3_lookupTable_[123].f90, with the module of shared
# procedures src/p3_lookupTable_common.f90) is compiled once (double-precision default reals,
# value-safe floating point) and run once per block of its
# outer loops; with 'all', the blocks of every table share one pool.  The blocks are assembled
# in index order (LT2: the four (Fl1,Fl2) runs are checked row by row and joined in the order
# p3_init reads).  Before a table is written, the script checks that every block exited 0,
# that all blocks of a table hold the same number of rows, that the header carries the
# generator's version and line 2 is empty, and that the row counts match the table's axes.
# The result does not depend on -j.
#
# Output (in --outdir): p3_lookupTable_1.dat-v<version>-2momI and -3momI,
# p3_lookupTable_2.dat-v<version>, p3_lookupTable_3.dat-v<version> (as requested; <version>
# is the generator's 'version' parameter), each with a .info file (code commit, compiler,
# flags, host, blocks, time).  A failed run writes no table (and removes an older one of the
# same name); its blocks' logs stay in DIR/work-<table>/log_*.  One run at a time per DIR.
#
# Cost (one core per block; measured with ifx on a 2025 x86-64 server, other machines differ):
#   blocks x time per block = core-hours
#   lt1-3mom 220 x ~30 min = ~115     lt1-2mom 20 x ~15 min = ~5
#   lt2      100 x ~7.5 min = ~12.5   lt3      400 x ~15 s  = ~1.5
#   all: ~135 core-hours, i.e. 40 min on 128 cores, ~14 h on 10 cores.
#
# Runs anywhere with bash (3.2 or later), a Fortran compiler and xargs (Linux, macOS).  On a
# cluster, run it inside a batch job on one node, with -j up to the cores requested.
set -euo pipefail

SRC=$(cd "$(dirname "$0")" && pwd)
usage() { awk 'NR > 1 && /^#/ {sub(/^# ?/, ""); print; next} NR > 1 {exit}' "$0"; exit "${1:-2}"; }   # the header above

TABLE=${1:-}; [ -n "$TABLE" ] || usage; case $TABLE in -h|--help) usage 0 ;; esac; shift
NJ=; FC=; OUT=$PWD; KEEP=0
while [ $# -gt 0 ]; do
  case $1 in
    -j) NJ=$2; shift 2 ;;
    --fc) FC=$2; shift 2 ;;
    --outdir) OUT=$2; shift 2 ;;
    --keep) KEEP=1; shift ;;
    -h|--help) usage 0 ;;
    *) echo "make_luts.sh: unknown option $1" >&2; usage ;;
  esac
done
case $TABLE in lt1-2mom|lt1-3mom|lt2|lt3|all) ;; *) echo "make_luts.sh: unknown table '$TABLE'" >&2; usage ;; esac
[ -n "$NJ" ] || NJ=1
if [ -z "$FC" ]; then
  for c in ifx ifort gfortran; do command -v $c >/dev/null 2>&1 && { FC=$c; break; }; done
  [ -n "$FC" ] || { echo "make_luts.sh: no Fortran compiler found (ifx, ifort, gfortran)" >&2; exit 1; }
fi
case $(basename "$FC") in
  ifx*|ifort*) FFLAGS="-r8 -O2 -fp-model precise" ;;
  gfortran*)   FFLAGS="-fdefault-real-8 -fdefault-double-8 -O2 -ffp-contract=off" ;;
  *) echo "make_luts.sh: unsupported compiler '$FC'" >&2; exit 1 ;;
esac
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd); export OUT
# one run per output directory at a time (mkdir is atomic)
mkdir "$OUT/.make_luts.lock" 2>/dev/null \
  || { echo "make_luts.sh: another run is using $OUT (remove $OUT/.make_luts.lock if none is)" >&2; exit 1; }
trap 'rmdir "$OUT/.make_luts.lock" 2>/dev/null || true' EXIT

log() { echo "[$(date -u +%H:%M:%S)] $*"; }
code_id() {  # commit of the generators, '-dirty' if src/ has uncommitted changes, and the branch
  git -C "$SRC" rev-parse --short HEAD >/dev/null 2>&1 || { echo 'not a git checkout'; return; }
  echo "$(git -C "$SRC" rev-parse --short HEAD)$(git -C "$SRC" diff --quiet HEAD -- . || echo -dirty) ($(git -C "$SRC" rev-parse --abbrev-ref HEAD))"
}

# --- per-table definitions: generator, table version (the header's version token), block
#     arguments (one block per line, in table order), block file name, final file name, and
#     SHAPE = the expected data rows as 'columns:rows' pairs (LT1 has a main table and an
#     ice-rain sub-table of a different width)
version_of() { sed -n "s/^ *character(len=[0-9]*), *parameter *:: *version *= *'\([^']*\)'.*/\1/p" "$SRC/$1" | head -1; }

def_lt1_3mom() {
  GEN=create_p3_lookupTable_1.f90; VER=$(version_of $GEN)-3momI
  FINAL=p3_lookupTable_1.dat-v$VER
  blocks() { for z in $(seq 1 11); do for r in 1 2 3 4 5; do for f in 1 2 3 4; do echo "3 $z $r $f"; done; done; done; }
  blockfile() { printf 'lt1-3mom_%02d_%d_%d.dat' "$2" "$3" "$4"; }
  SHAPE="36:$((11*5*4*4*50)) 7:$((11*5*4*4*50*30))"
}
def_lt1_2mom() {
  GEN=create_p3_lookupTable_1.f90; VER=$(version_of $GEN)-2momI
  FINAL=p3_lookupTable_1.dat-v$VER
  blocks() { for r in 1 2 3 4 5; do for f in 1 2 3 4; do echo "2 $r $f"; done; done; }
  blockfile() { printf 'lt1-2mom_%d_%d.dat' "$2" "$3"; }
  SHAPE="25:$((5*4*4*50)) 6:$((5*4*4*50*30))"
}
def_lt2() {
  GEN=create_p3_lookupTable_2.f90; VER=$(version_of $GEN)
  FINAL=p3_lookupTable_2.dat-v$VER
  blocks() { for fl in "0 0" "0 1" "1 0" "1 1"; do for i in $(seq 1 25); do echo "$fl $i"; done; done; }
  blockfile() { printf 'lt2_%d%d_%02d.dat' "$1" "$2" "$3"; }
  SHAPE="14:$((25*4*5*25*4*5))"
  assemble() { assemble_lt2 "$@"; }
}
def_lt3() {
  GEN=create_p3_lookupTable_3.f90; VER=$(version_of $GEN)
  FINAL=p3_lookupTable_3.dat-v$VER
  blocks() { for z in $(seq 1 80); do for r in 1 2 3 4 5; do echo "$z $r"; done; done; }
  blockfile() { printf 'lt3_%02d_%d.dat' "$1" "$2"; }
  SHAPE="9:$((80*5*4*4*50))"
}

# prepare_table T: compile the generator in a fresh work dir and append its blocks to the pool
prepare_table() {
  local T=$1 W nb
  assemble() { assemble_cat "$@"; }   # default; a table definition may replace it
  def_${T//-/_}
  W=$OUT/work-$T; rm -rf "$W"; mkdir -p "$W"
  rm -f "$OUT/$FINAL" "$OUT/$FINAL.info"     # never leave an older table beside a failed run
  ( cd "$W" && $FC $FFLAGS -o gen "$SRC/p3_lookupTable_common.f90" "$SRC/$GEN" ) > "$W/compile.log" 2>&1 \
    || { cat "$W/compile.log"; echo "make_luts.sh: $T: compile failed ($W/compile.log)" >&2; exit 1; }
  blocks > "$W/blocks.txt"; nb=$(wc -l < "$W/blocks.txt" | tr -d ' ')
  sed "s|^|$T |" "$W/blocks.txt" >> "$OUT/pool.txt"
  log "$T: $GEN (table version $VER), $nb blocks, work dir $W"
}

# stop_pool: stop xargs (no new blocks), then the running blocks (xargs -> sh -> gen)
stop_pool() {
  local kids p
  kids=$(pgrep -P "$XPID" 2>/dev/null || true)
  kill "$XPID" 2>/dev/null || true
  for p in $kids; do pkill -TERM -P "$p" 2>/dev/null || true; kill "$p" 2>/dev/null || true; done
  echo "make_luts.sh: stopped; no table written" >&2
}

# run_pool: run every block of every prepared table, -j at a time.  Each block runs in its
# table's work dir ($OUT/work-<table>), logs to its own file, and records its exit status.
run_pool() {
  local nb nfail
  nb=$(wc -l < "$OUT/pool.txt" | tr -d ' ')
  log "running $nb blocks, -j $NJ ($FC $FFLAGS)"
  T0=$(date +%s)
  # xargs runs in the background so that a TERM/INT sent to this script alone stops the blocks too
  xargs -P "$NJ" -L 1 sh -c 'cd "$OUT/work-$0" && ./gen "$@" > "log_$(echo "$@" | tr " " _)" 2>&1; echo "$? $0 $*"' \
    < "$OUT/pool.txt" > "$OUT/pool_status.txt" &
  XPID=$!
  trap 'stop_pool; exit 143' TERM INT
  wait "$XPID" || true
  trap - TERM INT
  T1=$(date +%s)
  nfail=$(awk '$1 != 0' "$OUT/pool_status.txt" | wc -l | tr -d ' ')
  log "blocks done in $((T1-T0)) s; failed: $nfail"
  if [ "$nfail" != 0 ] || [ "$(wc -l < "$OUT/pool_status.txt" | tr -d ' ')" != "$nb" ]; then
    awk '$1 != 0' "$OUT/pool_status.txt" > "$OUT/pool_failed.txt"
    sed -n '1,10p' "$OUT/pool_failed.txt"
    echo "make_luts.sh: block(s) failed or missing (status, table, block: $OUT/pool_failed.txt; logs in $OUT/work-*/log_*)" >&2
    exit 1
  fi
}

# finish_table T: assemble the blocks, check the table, write it and its provenance file
finish_table() {
  local T=$1 W nb
  assemble() { assemble_cat "$@"; }
  def_${T//-/_}
  W=$OUT/work-$T; nb=$(wc -l < "$W/blocks.txt" | tr -d ' ')
  check_blocks "$W"
  assemble "$W" "$W/table"
  check_table "$W/table"
  {
    echo "table:      $FINAL"
    echo "generator:  src/$GEN (table version $VER)"
    echo "code:       ${P3_LUT_CODE:-$(code_id)}"
    echo "compiler:   $($FC --version 2>&1 | head -1)"
    echo "flags:      $FFLAGS"
    echo "host:       $(hostname) ($(uname -sm))"
    echo "blocks:     $nb ($(wc -l < "$OUT/pool.txt" | tr -d ' ') in the run), -j $NJ, $((T1-T0)) s"
    echo "date:       $(date -u +%FT%TZ)"
  } > "$W/info"
  mv "$W/info" "$OUT/$FINAL.info"
  mv "$W/table" "$OUT/$FINAL"
  log "$T: wrote $OUT/$FINAL (+ .info)"
  [ "$KEEP" = 1 ] || rm -rf "$W"
}

# assemble WORKDIR OUT: concatenate the blocks in index order (the order of blocks.txt),
# never in directory-listing order
assemble_cat() {
  local W=$1 out=$2 a f
  : > "$out"
  while read -r a; do
    f=$(blockfile $a); [ -s "$W/$f" ] || { echo "make_luts.sh: missing $W/$f" >&2; exit 1; }
    cat "$W/$f" >> "$out"
  done < "$W/blocks.txt"
}

# check_blocks WORKDIR: every block of a table must hold the same rows (data rows counted by
# number of columns; the first block's 2 header lines are not data).  A block that stopped
# early, or wrote twice, is named here instead of surfacing only as a wrong total.
check_blocks() {
  local W=$1 a f ref
  : > "$W/block_shapes.txt"
  while read -r a; do
    f=$(blockfile $a); [ -s "$W/$f" ] || { echo "make_luts.sh: $T: missing block file $f" >&2; exit 1; }
    printf '%s|%s\n' "$(awk 'NF >= 5 {n[NF]++} END {for (c in n) print c ":" n[c]}' "$W/$f" | sort -n | tr '\n' ' ')" "$f" \
      >> "$W/block_shapes.txt"
  done < "$W/blocks.txt"
  ref=$(cut -d'|' -f1 "$W/block_shapes.txt" | sort | uniq -c | sort -rn | awk 'NR == 1 {sub(/^ *[0-9]+ /, ""); print}')
  awk -F'|' -v ref="$ref" '$1 != ref {print $2 ": " $1}' "$W/block_shapes.txt" > "$W/block_outliers.txt"
  if [ -s "$W/block_outliers.txt" ]; then
    sed -n '1,10p' "$W/block_outliers.txt"
    echo "make_luts.sh: $T: block(s) differ from the others (columns:rows $ref)" >&2; exit 1
  fi
}

# LT2: one file per (Fl1,Fl2) run, each with the 6 index columns (30 characters) and its
# 2 data columns (30 characters); check that the index fields agree on every row, then
# join the data columns in the order p3_init reads: (0,0),(0,1),(1,0),(1,1)
assemble_lt2() {
  local W=$1 out=$2 fl i f
  for fl in 00 01 10 11; do
    : > "$W/run_$fl.txt"
    for i in $(seq 1 25); do
      f=$(blockfile ${fl:0:1} ${fl:1:1} $i); [ -s "$W/$f" ] || { echo "make_luts.sh: missing $W/$f" >&2; exit 1; }
      cat "$W/$f" >> "$W/run_$fl.txt"
    done
  done
  head -2 "$W/run_00.txt" > "$out"
  tail -n +3 "$W/run_00.txt" > "$W/data_00.txt"
  paste -d '\t' "$W/data_00.txt" "$W/run_01.txt" "$W/run_10.txt" "$W/run_11.txt" | awk -F '\t' '
    { for (k = 1; k <= 4; k++) if (length($k) != 60) {
        printf "make_luts.sh: LT2 row %d: a (Fl1,Fl2) run row is not 60 characters (format 6i5,2e15.6)\n", NR > "/dev/stderr"; bad = 1; exit 1 }
      for (k = 2; k <= 4; k++) if (substr($k, 1, 30) != substr($1, 1, 30)) {
        printf "make_luts.sh: LT2 row %d: index columns of the (Fl1,Fl2) runs disagree\n", NR > "/dev/stderr"; bad = 1; exit 1 }
      print substr($1, 1, 60) substr($2, 31, 30) substr($3, 31, 30) substr($4, 31, 30) }
    END { if (bad) exit 1 }' >> "$out" || exit 1
}

check_table() {  # header version token, an empty second line, and the row count of each width
  local f=$1 h got
  h=$(head -1 "$f" | awk '{print $NF}')
  [ "$h" = "$VER" ] || { echo "make_luts.sh: header version '$h' is not '$VER'" >&2; exit 1; }
  [ -z "$(sed -n 2p "$f")" ] || { echo "make_luts.sh: line 2 is not empty" >&2; exit 1; }
  got=$(tail -n +3 "$f" | awk '{n[NF]++} END {for (c in n) print c ":" n[c]}' | sort -n | tr '\n' ' ' | sed 's/ $//')
  [ "$got" = "$(echo $SHAPE | tr ' ' '\n' | sort -n | tr '\n' ' ' | sed 's/ $//')" ] \
    || { echo "make_luts.sh: data rows (columns:rows) are '$got', expected '$SHAPE'" >&2; exit 1; }
  log "$T: checked header '$VER' and rows (columns:rows) $got"
}

# all tables share one pool of blocks, so 'all' takes as long as its slowest block
if [ "$TABLE" = all ]; then TABLES="lt1-2mom lt1-3mom lt2 lt3"; else TABLES=$TABLE; fi
: > "$OUT/pool.txt"
for t in $TABLES; do prepare_table $t; done
run_pool
for t in $TABLES; do finish_table $t; done
[ "$KEEP" = 1 ] || rm -f "$OUT/pool.txt" "$OUT/pool_status.txt" "$OUT/pool_failed.txt"
