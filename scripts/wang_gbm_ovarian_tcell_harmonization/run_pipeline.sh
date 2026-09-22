#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_pipeline.sh CONFIG_YAML RUN_DIR" >&2
  exit 2
fi

CONFIG_PATH="$1"
RUN_DIR="$2"
SCRIPT_DIR="$(cd "$(dirname "$CONFIG_PATH")" && pwd)"
SIF_PATH="$(awk '$1 == "sif_path:" {print $2}' "$CONFIG_PATH")"
REQUIRED_HOST="$(awk '$1 == "required_hostname:" {print $2}' "$CONFIG_PATH")"

if [[ "$(hostname -s)" != "$REQUIRED_HOST" ]]; then
  echo "This workflow must run directly on $REQUIRED_HOST; observed $(hostname -s)" >&2
  exit 3
fi
if [[ ! -f "$SIF_PATH" ]]; then
  echo "SIF not found: $SIF_PATH" >&2
  exit 4
fi

mkdir -p "$RUN_DIR"/{provenance,audit,annotated_objects,cell_metadata,sample_summaries,qc_tables,figures,association_results,report,logs}
date --iso-8601=seconds > "$RUN_DIR/provenance/pipeline_started_at.txt"
hostname -f > "$RUN_DIR/provenance/hostname.txt"
/usr/bin/apptainer --version > "$RUN_DIR/provenance/apptainer_version.txt"
find "$SCRIPT_DIR" -type f -not -name '.DS_Store' -print0 | sort -z | xargs -0 sha256sum > "$RUN_DIR/provenance/code_sha256.tsv"
{
  awk '/^[[:space:]]+(GBM|Ovarian):/ {print $2}' "$CONFIG_PATH"
  printf '%s\n' "$SIF_PATH"
} | xargs -n1 sha256sum > "$RUN_DIR/provenance/input_and_sif_sha256.tsv"
cp "$CONFIG_PATH" "$RUN_DIR/provenance/config.yaml"
cp "$SCRIPT_DIR/marker_programs.tsv" "$RUN_DIR/provenance/marker_programs.tsv"

run_r() {
  local script="$1"
  echo "[$(date --iso-8601=seconds)] START $script"
  /usr/bin/apptainer exec --cleanenv \
    --bind /share/lab_crd:/share/lab_crd \
    --env OMP_NUM_THREADS=8 \
    --env OPENBLAS_NUM_THREADS=8 \
    --env MKL_NUM_THREADS=8 \
    "$SIF_PATH" Rscript "$SCRIPT_DIR/$script" "$CONFIG_PATH" "$RUN_DIR"
  echo "[$(date --iso-8601=seconds)] END $script"
}

run_r 02_harmonize_annotations.R
run_r 03_qc_and_comparison.R
run_r 04_fitness_association.R
run_r 05_build_report.R
date --iso-8601=seconds > "$RUN_DIR/provenance/pipeline_complete.txt"
echo "Pipeline completed successfully: $RUN_DIR"
