#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "${1:-$(dirname "${BASH_SOURCE[0]}")/..}" && pwd)"
ENV_PREFIX="${ENV_PREFIX:-${BASE_DIR}/env/alignment_crispresso}"
MAMBA_BIN="${MAMBA_BIN:-micromamba}"
N_PROCESSES="${N_PROCESSES:-8}"
cd "${BASE_DIR}"

if ! command -v "${MAMBA_BIN}" >/dev/null 2>&1; then
  echo "[ERROR] micromamba not found: ${MAMBA_BIN}" >&2
  exit 1
fi

LOG_DIR="${BASE_DIR}/logs"
QC_DIR="${BASE_DIR}/results/qc"
CRISPRESSO_DIR="${BASE_DIR}/results/crispresso"
PROCESSED_DIR="${BASE_DIR}/results/processed_forward_reads"

mkdir -p "${LOG_DIR}" "${QC_DIR}/fastqc" "${QC_DIR}/multiqc" "${CRISPRESSO_DIR}" "${PROCESSED_DIR}"

timestamp="$(date +%Y%m%d_%H%M%S)"
log_file="${LOG_DIR}/run_crispresso_workflow_${timestamp}.log"
exec > >(tee -a "${log_file}") 2>&1

echo "[INFO] Start: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "[INFO] Base directory: ${BASE_DIR}"
echo "[INFO] Environment prefix: ${ENV_PREFIX}"
echo "[INFO] N processes: ${N_PROCESSES}"

FASTQ_DIR="${BASE_DIR}/data/fastq"
if [[ ! -d "${FASTQ_DIR}" ]]; then
  echo "[ERROR] Raw FASTQ directory missing: ${FASTQ_DIR}" >&2
  exit 1
fi
if [[ ! -d "${ENV_PREFIX}" ]]; then
  echo "[ERROR] CRISPResso environment missing: ${ENV_PREFIX}" >&2
  exit 1
fi

echo "[INFO] Extracting forward-oriented informative reads"
while IFS=$'\t' read -r specimen_id cohort locus sample_name r1 r2 remote_r1 remote_r2 analysis_read analysis_fastq analysis_name processed_summary_json fastq_status analysis_include exclusion_reason; do
  if [[ "${analysis_name}" == "analysis_name" ]]; then
    continue
  fi
  if [[ "${analysis_include}" != "yes" ]]; then
    echo "[INFO] Skipping excluded library ${analysis_name} (${exclusion_reason})"
    continue
  fi

  for input_fastq in "${remote_r1}" "${remote_r2}"; do
    if [[ ! -f "${input_fastq}" ]]; then
      echo "[ERROR] Required FASTQ missing for ${analysis_name}: ${input_fastq}" >&2
      exit 1
    fi
  done

  primer_fwd="$(awk -F'\t' -v locus="${locus}" 'NR>1 && $1==locus {print $4}' "${BASE_DIR}/references/amplicons.tsv")"
  "${MAMBA_BIN}" run -p "${ENV_PREFIX}" \
    python3 "${BASE_DIR}/source/extract_forward_informative_reads.py" \
    --fastq-r1 "${remote_r1}" \
    --fastq-r2 "${remote_r2}" \
    --primer-fwd "${primer_fwd}" \
    --output-fastq "${analysis_fastq}" \
    --summary-json "${processed_summary_json}"
done < "${BASE_DIR}/metadata/analysis_manifest.tsv"

echo "[INFO] Running FastQC on available FASTQ files"
"${MAMBA_BIN}" run -p "${ENV_PREFIX}" \
  fastqc \
  -t "${N_PROCESSES}" \
  -o "${QC_DIR}/fastqc" \
  "${FASTQ_DIR}"/*.fastq.gz \
  "${PROCESSED_DIR}"/*.fastq.gz

echo "[INFO] Running MultiQC"
"${MAMBA_BIN}" run -p "${ENV_PREFIX}" \
  multiqc \
  "${QC_DIR}/fastqc" \
  -o "${QC_DIR}/multiqc" \
  -n "alignment_multiqc_report.html"

echo "[INFO] Running CRISPRessoBatch for Kras (single-end R1; HDR-enabled)"
"${MAMBA_BIN}" run -p "${ENV_PREFIX}" \
  CRISPRessoBatch \
  --batch_settings "${BASE_DIR}/metadata/crispresso_kras_batch.tsv" \
  --batch_output_folder "${CRISPRESSO_DIR}/kras_batch" \
  --skip_failed \
  --place_report_in_output_folder \
  --write_detailed_allele_table \
  --suppress_amplicon_name_truncation \
  --exclude_bp_from_left 0 \
  --exclude_bp_from_right 0 \
  --quantification_window_size 0 \
  --quantification_window_center -3 \
  -p "${N_PROCESSES}"

echo "[INFO] Running CRISPRessoBatch for Trp53 (single-end R1; indel-focused)"
"${MAMBA_BIN}" run -p "${ENV_PREFIX}" \
  CRISPRessoBatch \
  --batch_settings "${BASE_DIR}/metadata/crispresso_trp53_batch.tsv" \
  --batch_output_folder "${CRISPRESSO_DIR}/trp53_batch" \
  --skip_failed \
  --place_report_in_output_folder \
  --write_detailed_allele_table \
  --suppress_amplicon_name_truncation \
  --exclude_bp_from_left 0 \
  --exclude_bp_from_right 0 \
  --quantification_window_size 5 \
  --quantification_window_center -3 \
  --ignore_substitutions \
  -p "${N_PROCESSES}"

echo "[INFO] Completed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
