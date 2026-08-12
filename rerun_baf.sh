#!/usr/bin/env bash
#SBATCH --job-name=rerun_baf
#SBATCH --output=slurm-rerun_baf-%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#
# Standalone rerun of just the pseudo-BAF step from run_cnv_ploidy_pipeline.sh.
# The CNVkit steps (coverage/segmentation/call) already completed successfully
# for all three samples -- no need to redo those. This only reruns section 4
# (pseudo-BAF at common SNP sites), using the CONTIG-FILTERED SNPs VCF to fix
# the earlier "sequence chr1_KI270766v1_alt not found" error (COMMON_SNPS_VCF
# had ALT-contig positions that hg38.analysisSet.fa, a no-alt reference,
# doesn't have).

set -euo pipefail

module load BCFtools
module load SAMtools

PROJECT_DIR="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/GC.PS.1929.WGS"
BAM_DIR="${PROJECT_DIR}/bams"
OUT_DIR="${PROJECT_DIR}/analysis/cnv_ploidy"
REF_FASTA="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/refdata/hg38.analysisSet.fa"

# the contig-filtered version -- adjust path if you saved it elsewhere
COMMON_SNPS_VCF="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/GC.PS.1929.WGS/vcf/gnomad.common_biallelic_snps.hg38.filtered.vcf.gz"

declare -A SAMPLES=(
  [parental]="rmdup_parental_deepWGS.bam"
  [resistant]="rmdup_resistant_deepWGS.bam"
  [fused]="rmdup_fused_deepWGS.bam"
)

mkdir -p "${OUT_DIR}/baf"

for name in "${!SAMPLES[@]}"; do
  bam="${BAM_DIR}/${SAMPLES[$name]}"
  echo "[baf] ${name}: ${bam}"

  bcftools mpileup -f "${REF_FASTA}" -R "${COMMON_SNPS_VCF}" -a AD -Ou "${bam}" \
    | bcftools call -m -Oz -o "${OUT_DIR}/baf/${name}.vcf.gz"

  tabix -p vcf "${OUT_DIR}/baf/${name}.vcf.gz"

  echo "[baf] ${name}: done, $(bcftools view -H "${OUT_DIR}/baf/${name}.vcf.gz" | wc -l) sites"
done

echo "Done. BAF vcfs in ${OUT_DIR}/baf/ -- rerun notebooks/cnv_ploidy_postprocessing.ipynb section 4."
