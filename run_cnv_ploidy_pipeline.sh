#!/usr/bin/env bash
#SBATCH --job-name=cnv_ploidy
#SBATCH --output=slurm-cnv_ploidy-%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#
# CNV + ploidy pipeline for GC.PS.1929.WGS (Parental / Resistant / Fused NSCLC cell line WGS)
#
# Assumes you are running from the project root (GC.PS.1929.WGS/) where
# bams/, Fastq/, genome.bed, genome.target.bed, genome.antitarget.bed and
# reference.cnn already exist from the earlier alignment + CNVkit setup.
#
# Tumor-only, no matched normal available. Strategy:
#   1. Absolute CN/segmentation per sample against the existing flat CNVkit
#      reference (or an NA12878/GIAB bam if you set NORMAL_BAM below).
#   2. Purity fixed at 1.0 (pure cell line, no stromal admixture) -- only
#      ploidy is estimated.
#   3. Pseudo-BAF from common population SNP sites, called directly on each
#      tumor bam (no matched normal needed) -- for LOH / allele-specific view.
#   4. ichorCNA as an independent ploidy/purity check.
#   5. Relative CN of Resistant and Fused against Parental (the real
#      biological comparison for an isogenic series).
#
# Requires on PATH: cnvkit.py, samtools, bcftools, bgzip/tabix, Rscript (ichorCNA),
# readCounter (HMMcopy utils; only needed for the ichorCNA step).

set -euo pipefail

### ---- User-editable paths -------------------------------------------------

PROJECT_DIR="$(pwd)"
BAM_DIR="${PROJECT_DIR}/bams"
OUT_DIR="${PROJECT_DIR}/analysis/cnv_ploidy"
REF_FASTA="/path/to/GRCh38.fa"                 # same reference used for alignment
REFFLAT="/path/to/refFlat.hg38.txt"             # for cnvkit annotate, optional but recommended
COMMON_SNPS_VCF="/path/to/gnomad.common_biallelic_snps.hg38.vcf.gz"  # for pseudo-BAF

# Set this to a public diploid WGS bam (e.g. NA12878 GIAB, aligned identically
# to your samples) if you want a real (non-flat) CNVkit reference. Leave empty
# to keep using the flat reference.cnn that's already in PROJECT_DIR.
NORMAL_BAM=""

# ichorCNA panel-of-normals / support files (ships with the ichorCNA R package)
ICHOR_DIR="/path/to/ichorCNA"
ICHOR_GC_WIG="${ICHOR_DIR}/inst/extdata/gc_hg38_1000kb.wig"
ICHOR_MAP_WIG="${ICHOR_DIR}/inst/extdata/map_hg38_1000kb.wig"
ICHOR_CENTROMERE="${ICHOR_DIR}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt"
ICHOR_PON=""   # leave empty -- no PoN available; ichorCNA will run without one

declare -A SAMPLES=(
  [parental]="rmdup_parental_deepWGS.bam"
  [resistant]="rmdup_resistant_deepWGS.bam"
  [fused]="rmdup_fused_deepWGS.bam"
)

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUT_DIR}"/{cnvkit,baf,ichorcna,seg}

### ---- 1. CNVkit reference --------------------------------------------------

if [[ -n "${NORMAL_BAM}" ]]; then
  echo "[ref] building CNVkit reference from ${NORMAL_BAM}"
  cnvkit.py batch \
    -n "${NORMAL_BAM}" \
    -f "${REF_FASTA}" \
    -t "${PROJECT_DIR}/genome.target.bed" \
    -g "${PROJECT_DIR}/genome.antitarget.bed" \
    --annotate "${REFFLAT}" \
    --output-reference "${OUT_DIR}/cnvkit/reference.cnn" \
    -p "${THREADS}" \
    --output-dir "${OUT_DIR}/cnvkit"
  CNVKIT_REF="${OUT_DIR}/cnvkit/reference.cnn"
else
  echo "[ref] reusing existing flat reference: ${PROJECT_DIR}/reference.cnn"
  CNVKIT_REF="${PROJECT_DIR}/reference.cnn"
fi

### ---- 2. Per-sample coverage, segmentation, calling ------------------------

for name in "${!SAMPLES[@]}"; do
  bam="${BAM_DIR}/${SAMPLES[$name]}"
  echo "[cnvkit] ${name}: ${bam}"

  cnvkit.py batch "${bam}" \
    -r "${CNVKIT_REF}" \
    -p "${THREADS}" \
    -d "${OUT_DIR}/cnvkit" \
    --scatter --diagram

  cnr="${OUT_DIR}/cnvkit/$(basename "${bam}" .bam).cnr"
  cns="${OUT_DIR}/cnvkit/$(basename "${bam}" .bam).cns"

  # Absolute copy number: purity fixed at 1.0 (pure cell line).
  # Ploidy left at default 2 here -- the notebook re-derives an empirical
  # ploidy estimate from the segment log2 distribution and you can re-run
  # this `call` step with --ploidy <N> once you have that number (or the
  # flow cytometry DNA-index result).
  cnvkit.py call "${cns}" \
    --purity 1.0 \
    --ploidy 2 \
    -o "${OUT_DIR}/cnvkit/${name}.call.cns"

  # SEG file for IGV / notebook
  cnvkit.py export seg "${cns}" -o "${OUT_DIR}/seg/${name}.seg"

  cp "${cnr}" "${OUT_DIR}/cnvkit/${name}.cnr"
  cp "${cns}" "${OUT_DIR}/cnvkit/${name}.cns"
done

### ---- 3. Relative CN vs. Parental (isogenic comparison) -------------------
# CNVkit can directly diff two .cnr files' log2 ratios by using one sample's
# .cnn-like coverage as the "reference" for the other. Simpler + more
# transparent here: just diff the log2 columns per-bin in the notebook using
# the .cnr files produced above (parental.cnr, resistant.cnr, fused.cnr).
# No extra command needed at this stage -- see notebook section 4.

### ---- 4. Pseudo-BAF at common SNP sites (no matched normal needed) --------

for name in "${!SAMPLES[@]}"; do
  bam="${BAM_DIR}/${SAMPLES[$name]}"
  echo "[baf] ${name}"
  bcftools mpileup -f "${REF_FASTA}" -R "${COMMON_SNPS_VCF}" -a AD -Ou "${bam}" \
    | bcftools call -m -Oz -o "${OUT_DIR}/baf/${name}.vcf.gz"
  tabix -p vcf "${OUT_DIR}/baf/${name}.vcf.gz"
done

### ---- 5. ichorCNA independent ploidy/purity estimate -----------------------

for name in "${!SAMPLES[@]}"; do
  bam="${BAM_DIR}/${SAMPLES[$name]}"
  wig="${OUT_DIR}/ichorcna/${name}.wig"
  echo "[ichorCNA] ${name}: generating read counts"
  readCounter --window 1000000 --quality 20 \
    --chromosome "$(seq -s, 1 22),X" \
    "${bam}" > "${wig}"

  echo "[ichorCNA] ${name}: running"
  Rscript "${ICHOR_DIR}/scripts/runIchorCNA.R" \
    --id "${name}" \
    --WIG "${wig}" \
    --gcWig "${ICHOR_GC_WIG}" \
    --mapWig "${ICHOR_MAP_WIG}" \
    --centromere "${ICHOR_CENTROMERE}" \
    ${ICHOR_PON:+--normalPanel "${ICHOR_PON}"} \
    --normal "0" \
    --scStates "c(1,3)" \
    --ploidy "c(2,3,4)" \
    --maxCN 8 \
    --outDir "${OUT_DIR}/ichorcna"
done

echo "Done. Outputs in ${OUT_DIR}/{cnvkit,baf,ichorcna,seg}"
echo "Next: open notebooks/cnv_ploidy_postprocessing.ipynb, point INPUT_DIR at ${OUT_DIR}"
