#!/usr/bin/env bash
#SBATCH --job-name=mark_dup
#SBATCH --output=slurm-mark_dup-%j.out
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#
# Duplicate-marking for GC.PS.1929.WGS bams (Parental / Resistant / Fused).
# Input bams (fused_deepWGS.bam, parental_deepWGS.bam, resistant_deepWGS.bam)
# came straight out of `bwa mem | samtools sort` with no dedup step, so they
# need this before any coverage-based CNV calling -- PCR/optical duplicates
# inflate coverage non-uniformly across bins and distort log2 ratios.
#
# Uses samtools markdup (needs name-collated + fixmate'd input first -- the
# standard samtools workflow), chained through pipes so no full-size
# intermediate bam ever hits disk. Produces rmdup_<name>_deepWGS.bam,
# matching the naming already used elsewhere in this project.
#
# module load SAMtools   # match whatever module gave you `samtools` for alignment

set -euo pipefail

BAM_DIR="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/GC.PS.1929.WGS/bams"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

declare -A SAMPLES=(
  [parental]="parental_deepWGS.bam"
  [resistant]="resistant_deepWGS.bam"
  [fused]="fused_deepWGS.bam"
)

for name in "${!SAMPLES[@]}"; do
  in_bam="${BAM_DIR}/${SAMPLES[$name]}"
  out_bam="${BAM_DIR}/rmdup_${SAMPLES[$name]}"

  echo "[markdup] ${name}: ${in_bam} -> ${out_bam}"

  samtools collate -@ "${THREADS}" -O -u "${in_bam}" \
    | samtools fixmate -@ "${THREADS}" -m -u - - \
    | samtools sort -@ "${THREADS}" -u - \
    | samtools markdup -@ "${THREADS}" -r -s - "${out_bam}"

  samtools index -@ "${THREADS}" "${out_bam}"

  echo "[markdup] ${name}: done, $(samtools view -c "${out_bam}") reads remaining"
done

echo "Done. rmdup_*.bam files are in ${BAM_DIR}, ready for run_cnv_ploidy_pipeline.sh"
