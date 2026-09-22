#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_direct_hpctpa3pc0028.sh RUN_ID CODE_DIR" >&2
  exit 2
fi

RUN_ID="$1"
CODE_DIR="$2"
OUTPUT_ROOT="/share/lab_crd/taoli/Project/CarTKFLs/results/wang_gbm_ovarian_tcell_harmonization"
RUN_DIR="$OUTPUT_ROOT/$RUN_ID"
CONFIG_PATH="$CODE_DIR/config.yaml"

if [[ "$(hostname -s)" != "hpctpa3pc0028" ]]; then
  echo "Run this launcher directly on hpctpa3pc0028; observed $(hostname -s)" >&2
  exit 3
fi
if [[ ! -f "$RUN_DIR/runtime_library/UCell/DESCRIPTION" ]]; then
  echo "Run-specific UCell library is missing: $RUN_DIR/runtime_library/UCell" >&2
  exit 4
fi

mkdir -p "$RUN_DIR/logs" "$RUN_DIR/provenance"
ATTEMPT_TAG="$(date +%Y%m%d_%H%M%S)"
LOG_PATH="$RUN_DIR/logs/pipeline_attempt_${ATTEMPT_TAG}.log"
nohup bash "$CODE_DIR/run_pipeline.sh" "$CONFIG_PATH" "$RUN_DIR" \
  > "$LOG_PATH" 2>&1 < /dev/null &
PID=$!
printf '%s\n' "$PID" > "$RUN_DIR/provenance/pipeline.pid"
printf '%s\n' "$LOG_PATH" > "$RUN_DIR/provenance/current_pipeline_log.txt"
printf 'pid=%s\nrun_dir=%s\nlog=%s\n' "$PID" "$RUN_DIR" "$LOG_PATH"
