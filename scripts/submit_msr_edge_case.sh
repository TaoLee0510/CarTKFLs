#!/usr/bin/env bash
set -euo pipefail

PROJECT=/share/lab_crd/taoli/Project/CarTKFLs
PAN_REPO=/share/lab_crd/taoli/Project/PANcanKFLs
SIF_PATH=/share/lab_crd/taoli/Docker/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif
RUN_DIR=${PROJECT}/results/runs/20260917_f93b2d2_uncapped_f4633d2/correlation_only
CONFIG_PATH=${RUN_DIR}/downstream.yaml
QOS=xxlarge
TIME_LIMIT=12:00:00

if ! command -v sbatch >/dev/null 2>&1 || ! command -v apptainer >/dev/null 2>&1; then
  set +eu
  # shellcheck disable=SC1091
  source /etc/profile
  set -eu
fi
for name in ${!SBATCH_@}; do unset "${name}"; done
for path in "${CONFIG_PATH}" "${RUN_DIR}/manifests/msr_samples.tsv" \
            "${RUN_DIR}/analysis/results/MSR/CarT_high_cn_8/steadyStatePredictions.Rds"; do
  [[ -s ${path} ]] || { echo "Missing input: ${path}" >&2; exit 1; }
done
[[ ! -s ${RUN_DIR}/msr_edge_case.tsv ]] || {
  echo "MSR edge-case result already exists; refusing to resubmit" >&2
  exit 1
}
[[ $(find "${RUN_DIR}/status" -maxdepth 1 -name 'msr__*.tsv' -type f | wc -l) -eq 12 ]] || {
  echo "Expected 12 patient-level MSR statuses" >&2
  exit 1
}
[[ -s ${RUN_DIR}/status/msr_reduce__CarT_high_cn_5.tsv &&
   -s ${RUN_DIR}/status/msr_reduce__CarT_high_cn_6.tsv ]] || {
  echo "Completed high_cn_5 and high_cn_6 reductions are required" >&2
  exit 1
}

reduce_job=$(sbatch --parsable --job-name=CarTMSREdgeCase \
  --ntasks=1 --cpus-per-task=1 --mem=32G --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/msr_edge_case_%j.out" \
  --error="${RUN_DIR}/logs/msr_edge_case_%j.err" \
  --export="ALL,RUN_STAGE=msr_reduce_first_minimum,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
  "${PROJECT}/scripts/submit.sh")
reduce_job=${reduce_job%%;*}
printf '%s\t%s\t%s\t\t1\t1\t32G\t%s\t%s\t%s\t0\t\n' \
  "$(date --iso-8601=seconds)" msr_reduce_first_minimum "${reduce_job}" \
  "${TIME_LIMIT}" "${QOS}" msr_reduce_first_minimum >>"${RUN_DIR}/submissions.tsv"

audit_job=$(sbatch --parsable --dependency="afterany:${reduce_job}" \
  --job-name=CarTMSRAuditPearson --ntasks=1 --cpus-per-task=1 --mem=16G \
  --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/msr_audit_retry_%j.out" \
  --error="${RUN_DIR}/logs/msr_audit_retry_%j.err" \
  --export="ALL,RUN_STAGE=msr_audit,REPO_PATH=${PAN_REPO},CONFIG_PATH=${CONFIG_PATH},SIF_PATH=${SIF_PATH},RUN_DIR=${RUN_DIR},APPTAINER_BIND=${PROJECT}:${PROJECT},QOS=${QOS},TIME_LIMIT=${TIME_LIMIT}" \
  "${PAN_REPO}/scripts/hpc/hnsc_melanoma_ovarian_downstream/submit_workflow.sh")
audit_job=${audit_job%%;*}
printf '%s\t%s\t%s\t%s\t1\t1\t16G\t%s\t%s\t%s\t0\t\n' \
  "$(date --iso-8601=seconds)" msr_audit_retry "${audit_job}" \
  "afterany:${reduce_job}" "${TIME_LIMIT}" "${QOS}" msr_audit \
  >>"${RUN_DIR}/submissions.tsv"
echo "[msr-edge-case] reduce=${reduce_job} audit=${audit_job} run_dir=${RUN_DIR}"
