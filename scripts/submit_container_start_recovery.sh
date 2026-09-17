#!/usr/bin/env bash
# Submit only the verified pre-container failures from the uncapped CarT run.
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail

PROJECT=/share/lab_crd/taoli/Project/CarTKFLs
EXPECTED_RUN=${PROJECT}/results/runs/20260917_f93b2d2_uncapped_f4633d2
RUN_DIR=${1:?Usage: submit_container_start_recovery.sh RUN_DIR}
CONFIG_PATH=${PROJECT}/config/analysis.yaml
if [[ "${RUN_DIR}" != "${EXPECTED_RUN}" ]]; then
  echo "Recovery must use the original run directory: ${EXPECTED_RUN}" >&2
  exit 1
fi
command -v sbatch >/dev/null 2>&1 || { echo "sbatch unavailable" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "python3 unavailable" >&2; exit 1; }
for name in ${!SBATCH_@}; do unset "${name}"; done

python3 "${PROJECT}/scripts/build_container_start_recovery.py" "${RUN_DIR}" --write
cp --no-clobber "${RUN_DIR}/fit_audit.tsv" \
  "${RUN_DIR}/fit_audit.before_container_start_recovery.tsv"

for minobs in 5 10 20; do
  manifest="${RUN_DIR}/manifests/container_start_recovery_MINOBS_${minobs}.tsv"
  [[ -f "${manifest}" ]] || continue
  tasks=$(($(wc -l <"${manifest}") - 1))
  ((tasks > 0)) || { echo "Empty recovery manifest: ${manifest}" >&2; exit 1; }
  case "${minobs}" in
    5) mem=48G ;;
    10) mem=16G ;;
    20) mem=8G ;;
  esac
  job=$(sbatch --parsable --job-name="CarTStartRecover${minobs}" \
    --array="1-${tasks}" --ntasks=1 --cpus-per-task=1 --mem="${mem}" \
    --qos=xxlarge --time=12:00:00 \
    --output="${RUN_DIR}/logs/recover${minobs}_%A_%a.out" \
    --error="${RUN_DIR}/logs/recover${minobs}_%A_%a.err" \
    --export="ALL,RUN_STAGE=fit,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH},TASK_FILE=${manifest},TASK_OFFSET=0" \
    "${PROJECT}/scripts/submit.sh")
  job=${job%%;*}
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "container_start_recovery_MINOBS_${minobs}" \
    "${job}" "" "${tasks}" 1 "${mem}" "12:00:00" "xxlarge" \
    "fit" 0 "${manifest}" >>"${RUN_DIR}/submissions.tsv"
  echo "MINOBS_${minobs} recovery_job=${job} tasks=${tasks} mem=${mem}"
done
