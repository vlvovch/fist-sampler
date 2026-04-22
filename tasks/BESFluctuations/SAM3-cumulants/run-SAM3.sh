#!/usr/bin/env bash
# Run the BES-SAM3-cumulants binary for a given collision energy and ensemble.
#
# Usage:
#   ./run-SAM3.sh <energy> [ensemble] [centrality]
#     energy:     7.7 | 14.5 | 19.6 | 27 | 39 | 62.4 | 200
#     ensemble:   GCE | B | BQS        (default: GCE)
#     centrality: directory label       (default: C0-5)
#
# Runs indefinitely (nevents=-1); cumulants flush every 1000 events.
# Stop with scancel (on HPC) or Ctrl-C (locally).

set -euo pipefail

# cd to repo root (script lives at tasks/BESFluctuations/SAM3-cumulants/)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/../../.."

ENERGY="${1:?usage: $0 <energy> [ensemble: GCE|B|BQS] [centrality, default C0-5]}"
ENSEMBLE="${2:-GCE}"
CENTRALITY="${3:-C0-5}"

case "$ENSEMBLE" in
  GCE) FLAGS=(--Bcanonical=0 --Qcanonical=0 --Scanonical=0) ;;
  B)   FLAGS=(--Bcanonical=1 --Qcanonical=0 --Scanonical=0) ;;
  BQS) FLAGS=(--Bcanonical=1 --Qcanonical=1 --Scanonical=1) ;;
  *)   echo "ERROR: unknown ensemble '$ENSEMBLE'. Use GCE|B|BQS." >&2; exit 1 ;;
esac

BIN=build/tasks/BESFluctuations/SAM3-cumulants/BES-SAM3-cumulants
INPUT=tasks/BESFluctuations/input/input.AuAu.${ENERGY}.${CENTRALITY}.EVHRG
SURFACE=input/hydro/AuAu.${ENERGY}/${CENTRALITY}/surface_eps_0.26.dat
OUTDIR=results/SAM3
OUTBASE=${OUTDIR}/AuAu.${ENERGY}.${CENTRALITY}.EVHRG
LOG=${OUTBASE}.${ENSEMBLE}.log

[[ -x "$BIN"     ]] || { echo "ERROR: binary not built: $BIN"               >&2; exit 1; }
[[ -f "$INPUT"   ]] || { echo "ERROR: input file missing: $INPUT"          >&2; exit 1; }
[[ -f "$SURFACE" ]] || { echo "ERROR: hypersurface missing: $SURFACE"      >&2; exit 1; }
mkdir -p "$OUTDIR"

echo "AuAu @ sqrt(s_NN) = ${ENERGY} GeV  |  ensemble=${ENSEMBLE}  |  centrality=${CENTRALITY}"
echo "Binary:       $BIN"
echo "Input:        $INPUT"
echo "Hypersurface: $SURFACE"
echo "Output base:  $OUTBASE  (binary appends .${ENSEMBLE}.SAM3-{cumulants,corrected}.dat)"
echo "Log:          $LOG"
echo

exec "$BIN" \
  "$INPUT" \
  --hypersurface_file="$SURFACE" \
  --output_file="${OUTBASE}.dat" \
  "${FLAGS[@]}" \
  --nevents=-1 \
  > "$LOG" 2>&1
