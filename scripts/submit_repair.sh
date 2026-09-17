#!/usr/bin/env bash
set -euo pipefail

PROJECT=/share/lab_crd/taoli/Project/CarTKFLs
CONFIG_PATH=${PROJECT}/config/analysis.yaml
RUN_DIR=${1:?Usage: submit_repair.sh RUN_DIR}
QOS=xxlarge
TIME_LIMIT=12:00:00

if ! command -v sbatch >/dev/null 2>&1; then
  set +eu
  # shellcheck disable=SC1091
  source /etc/profile
  set -eu
fi
for name in ${!SBATCH_@}; do unset "${name}"; done
[[ -f "${RUN_DIR}/tasks_all.tsv" && -f "${RUN_DIR}/submissions.tsv" ]] || {
  echo "Not a prepared CarTKFLs run: ${RUN_DIR}" >&2; exit 1;
}

fit_jobs=$(awk -F '\t' '$2 ~ /^fit_MINOBS_/ {print $3}' "${RUN_DIR}/submissions.tsv")
[[ -n "${fit_jobs}" ]] || { echo "No fit array IDs in submission ledger" >&2; exit 1; }
for job in ${fit_jobs}; do
  if [[ -n "$(squeue -h -j "${job}")" ]]; then
    echo "Fit array ${job} is still queued/running; wait for terminal accounting" >&2
    exit 1
  fi
done
original_audit=$(awk -F '\t' '$2 == "audit" {print $3; exit}' "${RUN_DIR}/submissions.tsv")
if [[ -n "${original_audit}" && -n "$(squeue -h -j "${original_audit}")" ]]; then
  echo "Original audit ${original_audit} has not finished" >&2
  exit 1
fi

build_job=$(sbatch --wait --parsable --job-name=CarTRepairManifest \
  --ntasks=1 --cpus-per-task=1 --mem=4G --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/repair_manifest_%j.out" \
  --error="${RUN_DIR}/logs/repair_manifest_%j.err" \
  --export="ALL,RUN_STAGE=build_repair,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
  "${PROJECT}/scripts/submit.sh")
build_job=${build_job%%;*}
printf '%s\trepair_manifest\t%s\t\t1\t1\t4G\t%s\t%s\tbuild_repair\t0\t\n' \
  "$(date --iso-8601=seconds)" "${build_job}" "${TIME_LIMIT}" "${QOS}" \
  >>"${RUN_DIR}/submissions.tsv"
manifest="${RUN_DIR}/repair_tasks.tsv"
total=$(($(wc -l <"${manifest}") - 1))
[[ "${total}" -gt 0 ]] || { echo "No eligible xval-only repairs"; exit 0; }

max_array_size=$(scontrol show config 2>/dev/null | awk -F= \
  '/^[[:space:]]*MaxArraySize[[:space:]]*=/ && !seen {gsub(/[[:space:]]/, "", $2); print $2; seen=1}')
if [[ ! "${max_array_size}" =~ ^[0-9]+$ ]] || ((max_array_size <= 1)); then
  max_array_size=1001
fi
max_chunk=$((max_array_size - 1))
offset=0
repair_jobs=()
while ((offset < total)); do
  chunk=$((total - offset))
  ((chunk > max_chunk)) && chunk=${max_chunk}
  job=$(sbatch --parsable --job-name=CarTXvalRepair \
    --array="1-${chunk}%32" --ntasks=1 --cpus-per-task=1 --mem=8G \
    --qos="${QOS}" --time="${TIME_LIMIT}" \
    --output="${RUN_DIR}/logs/repair_%A_%a.out" \
    --error="${RUN_DIR}/logs/repair_%A_%a.err" \
    --export="ALL,RUN_STAGE=repair,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH},TASK_FILE=${manifest},TASK_OFFSET=${offset}" \
    "${PROJECT}/scripts/submit.sh")
  job=${job%%;*}
  repair_jobs+=("${job}")
  printf '%s\trepair\t%s\t\t%s\t1\t8G\t%s\t%s\trepair\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "${job}" "${chunk}" "${TIME_LIMIT}" "${QOS}" \
    "${offset}" "${manifest}" >>"${RUN_DIR}/submissions.tsv"
  offset=$((offset + chunk))
done

dependency="afterany:$(IFS=:; echo "${repair_jobs[*]}")"
audit_job=$(sbatch --parsable --dependency="${dependency}" \
  --job-name=CarTRepairAudit --ntasks=1 --cpus-per-task=1 --mem=16G \
  --qos="${QOS}" --time="${TIME_LIMIT}" \
  --output="${RUN_DIR}/logs/repair_audit_%j.out" \
  --error="${RUN_DIR}/logs/repair_audit_%j.err" \
  --export="ALL,RUN_STAGE=audit,RUN_DIR=${RUN_DIR},CONFIG_PATH=${CONFIG_PATH}" \
  "${PROJECT}/scripts/submit.sh")
audit_job=${audit_job%%;*}
printf '%s\trepair_audit\t%s\t%s\t1\t1\t16G\t%s\t%s\taudit\t0\t\n' \
  "$(date --iso-8601=seconds)" "${audit_job}" "${dependency}" \
  "${TIME_LIMIT}" "${QOS}" >>"${RUN_DIR}/submissions.tsv"
echo "[repair] run_dir=${RUN_DIR} tasks=${total} repair_jobs=${repair_jobs[*]} audit_job=${audit_job}"
