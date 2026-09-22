#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: check_background_run.sh RUN_DIR" >&2
  exit 2
fi
RUN_DIR="$1"
PID_FILE="$RUN_DIR/provenance/pipeline.pid"
if [[ ! -f "$PID_FILE" ]]; then
  echo "PID file missing: $PID_FILE" >&2
  exit 3
fi
PID="$(tr -d '[:space:]' < "$PID_FILE")"
LOG_PATH_FILE="$RUN_DIR/provenance/current_pipeline_log.txt"
if [[ -f "$LOG_PATH_FILE" ]]; then
  LOG_PATH="$(tr -d '\r\n' < "$LOG_PATH_FILE")"
else
  LOG_PATH="$RUN_DIR/logs/pipeline.log"
fi
if kill -0 "$PID" 2>/dev/null; then
  STATUS="RUNNING"
elif [[ -f "$RUN_DIR/provenance/pipeline_complete.txt" ]]; then
  STATUS="COMPLETE"
else
  STATUS="NOT_RUNNING_INCOMPLETE"
fi
printf 'status=%s\npid=%s\n' "$STATUS" "$PID"
printf 'log=%s\n' "$LOG_PATH"
tail -n 40 "$LOG_PATH" 2>/dev/null || true
