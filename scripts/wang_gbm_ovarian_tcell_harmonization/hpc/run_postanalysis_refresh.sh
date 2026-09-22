#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_postanalysis_refresh.sh CODE_DIR RUN_DIR" >&2
  exit 2
fi
CODE_DIR="$1"
RUN_DIR="$2"
SIF_PATH="$(awk '$1 == "sif_path:" {print $2}' "$CODE_DIR/config.yaml")"
if [[ "$(hostname -s)" != "hpctpa3pc0028" ]]; then
  echo "This refresh must run directly on hpctpa3pc0028" >&2
  exit 3
fi

run_r() {
  /usr/bin/apptainer exec --cleanenv \
    --bind /share/lab_crd:/share/lab_crd \
    "$SIF_PATH" Rscript "$CODE_DIR/$1" "$CODE_DIR/config.yaml" "$RUN_DIR"
}

run_r 04_fitness_association.R
run_r 05_build_report.R
date --iso-8601=seconds > "$RUN_DIR/provenance/postanalysis_refresh_complete.txt"
