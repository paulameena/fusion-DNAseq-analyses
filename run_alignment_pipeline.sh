#!/usr/bin/env bash
#SBATCH --job-name=wgs_align
#SBATCH --output=slurm-wgs_align-%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#
# WGS alignment pipeline for GC.PS.1929.WGS (Parental / Resistant / Fused NSCLC cell line WGS)
#
# Fastq -> bwa-mem2 -> coordinate-sorted, duplicate-marked, indexed BAM.
# Run from the project root -- this is the step that produces bams/*.bam
# consumed by run_cnv_ploidy_pipeline.sh (which expects, e.g.,
# bams/rmdup_parental_deepWGS.bam).
#
# Layout expected:
#   Fastq/<prefix>*_R1*.f(ast)q.gz
#   Fastq/<prefix>*_R2*.f(ast)q.gz
# One sample may span multiple lanes -- every R1/R2 pair matching a sample's
# prefix is aligned separately (its own @RG) and then merged before dedup.
#
# Requires on PATH: bwa-mem2 (or bwa), samtools.

set -euo pipefail

### ---- User-editable paths -------------------------------------------------

PROJECT_DIR="$(pwd)"
FASTQ_DIR="${PROJECT_DIR}/Fastq"
BAM_DIR="${PROJECT_DIR}/bams"
REF_FASTA="/path/to/GRCh38.fa"          # same reference the CNV pipeline uses
BWA_BIN="bwa-mem2"                       # or "bwa" if that's what's installed

# sample key -> fastq filename prefix (glob) and output bam basename suffix,
# matching the SAMPLES map in run_cnv_ploidy_pipeline.sh
declare -A SAMPLE_PREFIX=(
  [parental]="Parental"
  [resistant]="Resistant"
  [fused]="Fused"
)
declare -A SAMPLE_OUT=(
  [parental]="rmdup_parental_deepWGS.bam"
  [resistant]="rmdup_resistant_deepWGS.bam"
  [fused]="rmdup_fused_deepWGS.bam"
)

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${BAM_DIR}"

### ---- 0. Reference indices --------------------------------------------------

if [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "[ref] indexing ${REF_FASTA} with samtools faidx"
  samtools faidx "${REF_FASTA}"
fi

if [[ ! -f "${REF_FASTA}.bwt.2bit.64" && ! -f "${REF_FASTA}.bwt" ]]; then
  echo "[ref] building ${BWA_BIN} index for ${REF_FASTA}"
  "${BWA_BIN}" index "${REF_FASTA}"
fi

### ---- Per-sample alignment --------------------------------------------------

for name in "${!SAMPLE_PREFIX[@]}"; do
  prefix="${SAMPLE_PREFIX[$name]}"
  out_bam="${BAM_DIR}/${SAMPLE_OUT[$name]}"
  work_dir="${BAM_DIR}/${name}.lanes"
  mkdir -p "${work_dir}"

  echo "[align] ${name}: searching ${FASTQ_DIR} for prefix '${prefix}'"
  mapfile -t r1_files < <(find "${FASTQ_DIR}" -maxdepth 1 -type f \
    \( -name "${prefix}*_R1*.fastq.gz" -o -name "${prefix}*_R1*.fq.gz" \) | sort)

  if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "[align] ${name}: no R1 fastqs found matching '${prefix}*_R1*', skipping" >&2
    continue
  fi

  lane_bams=()
  lane_idx=0
  for r1 in "${r1_files[@]}"; do
    r2="${r1/_R1/_R2}"
    if [[ ! -f "${r2}" ]]; then
      echo "[align] ${name}: no mate found for ${r1} (expected ${r2}), skipping lane" >&2
      continue
    fi

    lane_idx=$((lane_idx + 1))
    lane_tag="lane${lane_idx}"

    # Pull instrument/run/flowcell/lane out of the first read's header
    # (illumina casava format: @<instrument>:<run>:<flowcell>:<lane>:...)
    header_fields="$(zcat -f "${r1}" | head -n1 | tr -d '@' | cut -d' ' -f1)"
    flowcell="$(cut -d: -f3 <<<"${header_fields}")"
    lane_num="$(cut -d: -f4 <<<"${header_fields}")"
    rg_id="${flowcell:-${name}}.${lane_num:-${lane_idx}}"
    rg_pu="${rg_id}"

    lane_bam="${work_dir}/${lane_tag}.sorted.bam"
    echo "[align] ${name} ${lane_tag}: ${r1} + ${r2} -> ${lane_bam}"

    "${BWA_BIN}" mem -t "${THREADS}" \
      -R "@RG\tID:${rg_id}\tSM:${name}\tLB:${name}\tPL:ILLUMINA\tPU:${rg_pu}" \
      "${REF_FASTA}" "${r1}" "${r2}" \
      | samtools sort -@ "${THREADS}" -o "${lane_bam}" -
    samtools index "${lane_bam}"

    lane_bams+=("${lane_bam}")
  done

  if [[ ${#lane_bams[@]} -eq 0 ]]; then
    echo "[align] ${name}: no usable lanes, skipping sample" >&2
    continue
  fi

  merged_bam="${work_dir}/merged.sorted.bam"
  if [[ ${#lane_bams[@]} -gt 1 ]]; then
    echo "[align] ${name}: merging ${#lane_bams[@]} lane bams"
    samtools merge -f -@ "${THREADS}" "${merged_bam}" "${lane_bams[@]}"
  else
    cp "${lane_bams[0]}" "${merged_bam}"
  fi

  echo "[dedup] ${name}: marking duplicates -> ${out_bam}"
  samtools sort -@ "${THREADS}" -n -o "${work_dir}/nsorted.bam" "${merged_bam}"
  samtools fixmate -@ "${THREADS}" -m "${work_dir}/nsorted.bam" "${work_dir}/fixmate.bam"
  samtools sort -@ "${THREADS}" -o "${work_dir}/csorted.bam" "${work_dir}/fixmate.bam"
  samtools markdup -@ "${THREADS}" "${work_dir}/csorted.bam" "${out_bam}"
  samtools index "${out_bam}"

  rm -rf "${work_dir}"
  echo "[align] ${name}: done -> ${out_bam}"
done

echo "Done. BAMs in ${BAM_DIR}/"
echo "Next: scripts/run_cnv_ploidy_pipeline.sh"
