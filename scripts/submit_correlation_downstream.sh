#!/usr/bin/env bash
set -euo pipefail

PROJECT=/share/lab_crd/taoli/Project/CarTKFLs
PAN_REPO=/share/lab_crd/taoli/Project/PANcanKFLs
SIF_PATH=/share/lab_crd/taoli/Docker/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif
CONFIG_PATH=${PROJECT}/config/analysis.yaml
SOURCE_RUN_DIR=${1:-${PROJECT}/results/runs/20260917_f93b2d2_uncapped_f4633d2}
RUN_DIR=${SOURCE_RUN_DIR}/correlation_only
QOS=xxlarge
TIME_LIMIT=12:00:00

if [[ ${SOURCE_RUN_DIR} != ${PROJECT}/results/runs/20260917_f93b2d2_uncapped_f4633d2 ]]; then
  echo "Unexpected source run: ${SOURCE_RUN_DIR}" >&2
  exit 1
fi
if ! command -v sbatch >/dev/null 2>&1 || ! command -v apptainer >/dev/null 2>&1; then
  set +eu
  # shellcheck disable=SC1091
  source /etc/profile
  set -eu
fi
command -v sbatch >/dev/null 2>&1 || { echo "sbatch unavailable" >&2; exit 1; }
command -v apptainer >/dev/null 2>&1 || { echo "apptainer unavailable" >&2; exit 1; }
for name in ${!SBATCH_@}; do unset "${name}"; done

for path in "${SOURCE_RUN_DIR}/fit_audit.tsv" \
            "${SOURCE_RUN_DIR}/manifests/selected_samples.tsv" \
            "${SOURCE_RUN_DIR}/analysis/results/OptimizedParameters/GSE296419_CarT_days/final_parameters_opt.Rds" \
            "${SIF_PATH}"; do
  [[ -s ${path} ]] || { echo "Missing input: ${path}" >&2; exit 1; }
done
[[ ! -e ${RUN_DIR}/manifests/selected_samples.tsv ]] || {
  echo "Correlation-only selection already exists: ${RUN_DIR}" >&2
  exit 1
}
mkdir -p "${RUN_DIR}/logs" "${RUN_DIR}/manifests" "${RUN_DIR}/status"
if [[ ! -e ${RUN_DIR}/submissions.tsv ]]; then
  printf 'submitted_at\tstage\tjob_id\tdependency\ttasks\tcpus_per_task\tmem\ttime\tqos\tworker_stage\ttask_offset\tmanifest\n' \
    >"${RUN_DIR}/submissions.tsv"
fi
record_submission() {
  local stage=$1 job_id=$2 dependency=$3 cpus=$4 mem=$5
  printf '%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\t%s\t%s\t0\t\n' \
    "$(date --iso-8601=seconds)" "${stage}" "${job_id}" "${dependency}" \
    "${cpus}" "${mem}" "${TIME_LIMIT}" "${QOS}" "${stage}" \
    >>"${RUN_DIR}/submissions.tsv"
}

select_job=$(sbatch --parsable --job-name=CarTSelectPearson \
  --ntasks=1 --cpus-per-task=2 --mem=16G \
  --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/select_pearson_%j.out" \
  --error="${RUN_DIR}/logs/select_pearson_%j.err" \
  --export="ALL,RUN_STAGE=select_correlation,RUN_DIR=${RUN_DIR},SOURCE_RUN_DIR=${SOURCE_RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
  "${PROJECT}/scripts/submit.sh")
select_job=${select_job%%;*}
record_submission select_correlation "${select_job}" "" 2 16G

downstream_job=$(sbatch --parsable --dependency="afterok:${select_job}" \
  --job-name=CarTDownstreamPearson --ntasks=1 --cpus-per-task=1 --mem=16G \
  --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/downstream_controller_%j.out" \
  --error="${RUN_DIR}/logs/downstream_controller_%j.err" \
  --export="ALL,RUN_STAGE=after_selection,REPO_PATH=${PAN_REPO},CONFIG_PATH=${RUN_DIR}/downstream.yaml,SIF_PATH=${SIF_PATH},RUN_DIR=${RUN_DIR},APPTAINER_BIND=${PROJECT}:${PROJECT},QOS=${QOS},TIME_LIMIT=${TIME_LIMIT}" \
  "${PAN_REPO}/scripts/hpc/hnsc_melanoma_ovarian_downstream/submit_workflow.sh")
downstream_job=${downstream_job%%;*}
record_submission downstream_controller "${downstream_job}" "afterok:${select_job}" 1 16G
echo "[correlation-only] source=${SOURCE_RUN_DIR} variant=${RUN_DIR} select=${select_job} downstream_controller=${downstream_job}"
