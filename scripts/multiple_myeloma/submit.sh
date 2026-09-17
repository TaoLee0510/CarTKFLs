#!/usr/bin/env bash
set -euo pipefail

PROJECT=/share/lab_crd/taoli/Project/CarTKFLs
PAN_REPO=/share/lab_crd/taoli/Project/PANcanKFLs
NUMBAT=/share/lab_crd/CarTData/Numbat/GSE210079
METADATA=/share/lab_crd/CarTData/GSE210079
SIF_PATH=/share/lab_crd/taoli/Docker/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif
CONFIG_PATH=${CONFIG_PATH:-${PROJECT}/config/gse210079_multiple_myeloma.yaml}
RUN_STAGE=${RUN_STAGE:-submit}
QOS=xxlarge
TIME_LIMIT=12:00:00

if ! command -v sbatch >/dev/null 2>&1 || ! command -v apptainer >/dev/null 2>&1; then
  set +eu
  # shellcheck disable=SC1091
  source /etc/profile
  set -eu
fi
command -v sbatch >/dev/null 2>&1 || { echo "sbatch unavailable" >&2; exit 1; }
command -v apptainer >/dev/null 2>&1 || { echo "apptainer unavailable" >&2; exit 1; }
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  for name in ${!SBATCH_@}; do unset "${name}"; done
fi

container_rscript() {
  apptainer exec --cleanenv --pwd "${PROJECT}" \
    --bind "${PROJECT}:${PROJECT}" \
    --bind "${PAN_REPO}:${PAN_REPO}:ro" \
    --bind "${NUMBAT}:${NUMBAT}:ro" \
    --bind "${METADATA}:${METADATA}:ro" \
    --bind "/share/lab_crd/taoli/Docker:/share/lab_crd/taoli/Docker:ro" \
    --env "OMP_NUM_THREADS=1" --env "OPENBLAS_NUM_THREADS=1" \
    --env "MKL_NUM_THREADS=1" --env "R_FUTURE_FORK_ENABLE=false" \
    --env "CART_RUN_DIR=${RUN_DIR:-}" \
    --env "SLURM_JOB_ID=${SLURM_JOB_ID:-}" \
    --env "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-}" \
    --env "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}" \
    "${SIF_PATH}" /usr/local/bin/pancankfls-rscript --vanilla "$@"
}

record_submission() {
  local stage=$1 id=$2 dependency=$3 tasks=$4 cpus=$5 mem=$6 manifest=${7:-}
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "${stage}" "${id}" "${dependency}" \
    "${tasks}" "${cpus}" "${mem}" "${TIME_LIMIT}" "${QOS}" \
    "${stage}" "0" "${manifest}" >>"${RUN_DIR}/submissions.tsv"
}

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  : "${RUN_DIR:?RUN_DIR required inside Slurm}"
  case "${RUN_STAGE}" in
    preflight)
      container_rscript "${PROJECT}/scripts/multiple_myeloma/preflight.R" \
        "${CONFIG_PATH}" "${RUN_DIR}"
      ;;
    fit)
      task_row=$((TASK_OFFSET + SLURM_ARRAY_TASK_ID))
      container_rscript "${PROJECT}/scripts/multiple_myeloma/fit_one.R" \
        "${CONFIG_PATH}" "${RUN_DIR}" "${TASK_FILE}" "${task_row}"
      ;;
    audit)
      container_rscript "${PROJECT}/scripts/multiple_myeloma/audit_fit.R" \
        "${CONFIG_PATH}" "${RUN_DIR}"
      manifest="${RUN_DIR}/manifests/kfl_samples.tsv"
      tasks=$(($(wc -l <"${manifest}") - 1))
      kfl_job=$(sbatch --parsable --job-name=MM_KFL \
        --array="1-${tasks}" --ntasks=1 --cpus-per-task=4 --mem=32G \
        --qos="${QOS}" --time="${TIME_LIMIT}" \
        --output="${RUN_DIR}/logs/kfl_%A_%a.out" \
        --error="${RUN_DIR}/logs/kfl_%A_%a.err" \
        --export="ALL,RUN_STAGE=kfl,RUN_DIR=${RUN_DIR},TASK_FILE=${manifest},CONFIG_PATH=${CONFIG_PATH}" \
        "${PROJECT}/scripts/multiple_myeloma/submit.sh")
      kfl_job=${kfl_job%%;*}
      record_submission kfl "${kfl_job}" "" "${tasks}" 4 32G "${manifest}"
      select_job=$(sbatch --parsable --dependency="afterok:${kfl_job}" \
        --job-name=MM_Select --ntasks=1 --cpus-per-task=2 --mem=32G \
        --qos="${QOS}" --time="${TIME_LIMIT}" \
        --output="${RUN_DIR}/logs/select_%j.out" \
        --error="${RUN_DIR}/logs/select_%j.err" \
        --export="ALL,RUN_STAGE=select,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
        "${PROJECT}/scripts/multiple_myeloma/submit.sh")
      select_job=${select_job%%;*}
      record_submission select "${select_job}" "afterok:${kfl_job}" 1 2 32G
      downstream_job=$(sbatch --parsable --dependency="afterok:${select_job}" \
        --job-name=MM_Downstream --ntasks=1 --cpus-per-task=1 --mem=16G \
        --qos="${QOS}" --time="${TIME_LIMIT}" \
        --output="${RUN_DIR}/logs/downstream_%j.out" \
        --error="${RUN_DIR}/logs/downstream_%j.err" \
        --export="ALL,RUN_STAGE=after_selection,REPO_PATH=${PAN_REPO},CONFIG_PATH=${RUN_DIR}/downstream.yaml,SIF_PATH=${SIF_PATH},RUN_DIR=${RUN_DIR},APPTAINER_BIND=${PROJECT}:${PROJECT},QOS=${QOS},TIME_LIMIT=${TIME_LIMIT}" \
        "${PAN_REPO}/scripts/hpc/hnsc_melanoma_ovarian_downstream/submit_workflow.sh")
      downstream_job=${downstream_job%%;*}
      record_submission downstream_controller "${downstream_job}" "afterok:${select_job}" 1 1 16G
      echo "[audit] kfl=${kfl_job} select=${select_job} downstream_controller=${downstream_job}"
      ;;
    kfl)
      container_rscript "${PAN_REPO}/scripts/hpc/hnsc_melanoma_ovarian_downstream/sample_worker.R" \
        kfl "${RUN_DIR}/downstream.yaml" "${RUN_DIR}" "${TASK_FILE}"
      ;;
    select)
      container_rscript "${PAN_REPO}/scripts/hpc/hnsc_melanoma_ovarian_downstream/reduce_worker.R" \
        select "${RUN_DIR}/downstream.yaml" "${RUN_DIR}"
      container_rscript "${PROJECT}/scripts/multiple_myeloma/summarize_selection.R" \
        "${CONFIG_PATH}" "${RUN_DIR}"
      ;;
    *) echo "Unknown RUN_STAGE=${RUN_STAGE}" >&2; exit 1 ;;
  esac
  exit 0
fi

[[ "${RUN_STAGE}" == submit ]] || { echo "RUN_STAGE must be submit outside Slurm" >&2; exit 1; }
[[ -d "${PROJECT}/.git" && -f "${CONFIG_PATH}" && -f "${SIF_PATH}" ]] || {
  echo "Project checkout, config, or SIF missing" >&2; exit 1;
}
RUN_ID=${RUN_ID:-$(date +%Y%m%d_%H%M%S)_$(git -C "${PROJECT}" rev-parse --short=8 HEAD)}
RUN_DIR="${PROJECT}/multiple_myeloma/runs/${RUN_ID}"
[[ ! -e "${RUN_DIR}" ]] || { echo "Run directory already exists: ${RUN_DIR}" >&2; exit 1; }
mkdir -p "${RUN_DIR}/logs" "${RUN_DIR}/manifests" "${RUN_DIR}/status"
printf 'submitted_at\tstage\tjob_id\tdependency\ttasks\tcpus_per_task\tmem\ttime\tqos\tworker_stage\ttask_offset\tmanifest\n' \
  >"${RUN_DIR}/submissions.tsv"

preflight_job=$(sbatch --wait --parsable --job-name=MM_Preflight \
  --ntasks=1 --cpus-per-task=1 --mem=16G --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/preflight_%j.out" \
  --error="${RUN_DIR}/logs/preflight_%j.err" \
  --export="ALL,RUN_STAGE=preflight,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
  "${PROJECT}/scripts/multiple_myeloma/submit.sh")
preflight_job=${preflight_job%%;*}
record_submission preflight "${preflight_job}" "" 1 1 16G
[[ -s "${RUN_DIR}/tasks_all.tsv" && -s "${RUN_DIR}/provenance.tsv" ]] || {
  echo "Preflight did not produce manifests; see ${RUN_DIR}/logs" >&2; exit 1;
}

max_array_size=$(scontrol show config 2>/dev/null | awk -F= \
  '/^[[:space:]]*MaxArraySize[[:space:]]*=/ && !seen {gsub(/[[:space:]]/, "", $2); print $2; seen=1}')
if [[ ! "${max_array_size}" =~ ^[0-9]+$ ]] || ((max_array_size <= 1)); then
  max_array_size=1001
fi
max_chunk=$((max_array_size - 1))
fit_jobs=()
for patient in P16 P32 P33; do
  for minobs in 5 10 20; do
    manifest="${RUN_DIR}/tasks_${patient}_MINOBS_${minobs}.tsv"
    total=$(($(wc -l <"${manifest}") - 1))
    ((total > 0)) || continue
    mem=$(awk -F '\t' -v p="${patient}" -v m="${minobs}" \
      'NR>1 && $1==p && $2==m {print $5}' "${RUN_DIR}/fit_resources.tsv")
    [[ "${mem}" =~ ^[0-9]+G$ ]] || { echo "Invalid memory tier: ${patient}/${minobs}" >&2; exit 1; }
    offset=0
    while ((offset < total)); do
      chunk=$((total - offset))
      ((chunk > max_chunk)) && chunk=${max_chunk}
      job=$(sbatch --parsable --job-name="MM_${patient}_M${minobs}" \
        --array="1-${chunk}" --ntasks=1 --cpus-per-task=1 --mem="${mem}" \
        --qos="${QOS}" --time="${TIME_LIMIT}" \
        --output="${RUN_DIR}/logs/fit_${patient}_m${minobs}_%A_%a.out" \
        --error="${RUN_DIR}/logs/fit_${patient}_m${minobs}_%A_%a.err" \
        --export="ALL,RUN_STAGE=fit,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH},TASK_FILE=${manifest},TASK_OFFSET=${offset}" \
        "${PROJECT}/scripts/multiple_myeloma/submit.sh")
      job=${job%%;*}
      fit_jobs+=("${job}")
      record_submission "fit_${patient}_MINOBS_${minobs}" "${job}" "" "${chunk}" 1 "${mem}" "${manifest}#offset=${offset}"
      offset=$((offset + chunk))
    done
  done
done
((${#fit_jobs[@]} > 0)) || { echo "No fit jobs submitted" >&2; exit 1; }
dependency="afterany:$(IFS=:; echo "${fit_jobs[*]}")"
audit_job=$(sbatch --parsable --dependency="${dependency}" \
  --job-name=MM_FitAudit --ntasks=1 --cpus-per-task=1 --mem=16G \
  --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/audit_%j.out" \
  --error="${RUN_DIR}/logs/audit_%j.err" \
  --export="ALL,RUN_STAGE=audit,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
  "${PROJECT}/scripts/multiple_myeloma/submit.sh")
audit_job=${audit_job%%;*}
record_submission audit "${audit_job}" "${dependency}" 1 1 16G
echo "[submit] run_dir=${RUN_DIR} fit_jobs=${fit_jobs[*]} audit_job=${audit_job} qos=${QOS} time=${TIME_LIMIT}"
