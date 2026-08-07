#!/usr/bin/env bash
#SBATCH --job-name=cnv_ploidy
#SBATCH --output=slurm-cnv_ploidy-%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=16:00:00
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

### ---- Cluster environment modules -----------------------------------------
# sbatch jobs start a fresh non-interactive shell -- it does NOT source your
# .bashrc, so any `module load` you rely on interactively will NOT carry into
# this job unless it's loaded explicitly here. If you see errors like
# "ImportError: libffi.so.8: cannot open shared object file" or similar
# missing-shared-library failures, it means Python/cnvkit/samtools/etc. are
# resolving to a toolchain module that isn't actually loaded in this job.
#
# Find and uncomment/adjust the modules that match what's on PATH for you
# interactively:
#   module avail Python 2>&1
#   module avail bcftools 2>&1
#   module avail samtools 2>&1
#
module load Python/3.11.5-GCCcore-13.2.0
module load BCFtools
module load SAMtools
module load R

### ---- User-editable paths -------------------------------------------------

PROJECT_DIR="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/GC.PS.1929.WGS"
BAM_DIR="${PROJECT_DIR}/bams"
OUT_DIR="${PROJECT_DIR}/analysis/cnv_ploidy"
REF_FASTA="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/refdata/hg38.analysisSet.fa"    # same reference used for alignment
REFFLAT=""    # for cnvkit annotate, optional but recommended
COMMON_SNPS_VCF="/mnt/vstor/SOM_CCCC_JGS25/shultesp/data/GC.PS.1929.WGS/vcf/gnomad.common_biallelic_snps.hg38.vcf.gz"  # for pseudo-BAF

# Set this to a public diploid WGS bam (e.g. NA12878 GIAB, aligned identically
# to your samples) if you want a real (non-flat) CNVkit reference. Leave empty
# to keep using the flat reference.cnn that's already in PROJECT_DIR.
NORMAL_BAM=""

# ichorCNA panel-of-normals / support files (ships with the ichorCNA R package)
ICHOR_DIR="~/tools/ichorCNA"
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

# cnvkit.py batch's per-bin coverage pass is the memory-heavy step (each
# worker process holds its own chunk of BAM-processing state) -- it has
# already OOM'd twice at 8 and 11 workers even after raising --mem. Keep its
# parallelism lower and independent of THREADS (which other, lighter steps
# below still use at full width) rather than raising --mem indefinitely.
CNVKIT_THREADS=4

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
    -p "${CNVKIT_THREADS}" \
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
    -p "${CNVKIT_THREADS}" \
    -d "${OUT_DIR}/cnvkit" \
    --scatter --diagram

  cnr="${OUT_DIR}/cnvkit/$(basename "${bam}" .bam).cnr"
  cns="${OUT_DIR}/cnvkit/$(basename "${bam}" .bam).cns"

  # `cnvkit.py batch` can die partway (e.g. a worker OOM-killed during the
  # per-bin coverage pass) WITHOUT returning a nonzero exit code or printing
  # any error -- it just silently produces no .cnr/.cns. Check explicitly
  # here and stop with a clear message, instead of limping into `call` with
  # a file that doesn't exist and getting a confusing pandas traceback.
  if [[ ! -s "${cnr}" || ! -s "${cns}" ]]; then
    echo "[FATAL] ${name}: cnvkit.py batch did not produce ${cnr} and/or ${cns}." >&2
    echo "        This step is memory-heavy (coverage across ~1.8M WGS bins)" >&2
    echo "        and most likely died from an OOM kill with no traceback." >&2
    echo "        Check 'sacct -j \$SLURM_JOB_ID --format=MaxRSS' against --mem," >&2
    echo "        and consider lowering -p/THREADS to reduce peak concurrent memory." >&2
    exit 1
  fi

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
