#!/usr/bin/env bash
set -euo pipefail

# 04_qc_splitted.sh — Run sample QC on each population-separated file for imputation upload

# Inputs:
# --pfile PATH         (required) base dataset prefix (from 01/02, hg38, prepped)
# --keep-samples  PATH (required) PATH to population-separated sample list (*.list; IID-only or "FID IID")
# --fa PATH            (optional) hg38 FASTA (bgzip+.fai) for sex-check
# --outdir DIR         (default qc_splitted_out)
# --snps-only yes|no   (default yes) whether to restrict to SNPs only
# --drop-ambig yes|no  (default no) whether to drop ambiguous A/T and C/G SNPs
# --hwe-command STR    (default empty) extra HWE filter to pass to plink2 (e.g. "1e-6 midp")
# --geno-thres FLOAT   (default 0.02) extra per-variant missingness filter to pass to plink2 (e.g. 0.02)
# --maf-thres FLOAT    (default 0.01) extra MAF filter to pass to plink2 (e.g. 0.01)
# --vcf-out yes|no (default no) whether to output final QCed data in VCF format
# --threads N          (default 2)

PFILE=""
KEEPSAMPLES=""
FA=""
OUTDIR="qc_splitted_out"
SNPS_ONLY="yes"
DROP_AMBIG="no"
HWE_COMMAND=""
GENO_THRES="0.02"
MAF_THRES="0.01"
THREADS=2
VCF_OUT="no"
VCF_CHR_PREFIX="chr"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfile) PFILE="$2"; shift 2 ;;
    --keep-samples) KEEPSAMPLES="$2"; shift 2 ;;
    --fa) FA="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --snps-only) SNPS_ONLY="$2"; shift 2 ;;
    --drop-ambig) DROP_AMBIG="$2"; shift 2 ;;
    --hwe-command) HWE_COMMAND="$2"; shift 2 ;;
    --geno-thres) GENO_THRES="$2"; shift 2 ;;
    --maf-thres) MAF_THRES="$2"; shift 2 ;;
    --vcf-out) VCF_OUT="$2"; shift 2 ;;
    --vcf-chr-prefix) VCF_CHR_PREFIX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) sed -n '1,40p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -n "$PFILE" ]] || { echo "ERROR: --pfile required"; exit 1; }
[[ -n "$KEEPSAMPLES" ]] || { echo "ERROR: --keep-samples required"; exit 1; }
[[ -f "$KEEPSAMPLES" ]] || { echo "ERROR: keep-samples file not found: $KEEPSAMPLES"; exit 1; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/qc_splitted.log"; exec > >(tee -a "$LOG") 2>&1
trap 'echo "[ERROR] Failed at line $LINENO" >&2' ERR

POP=$(basename "$KEEPSAMPLES" .list)
echo "[INFO] QC on population-separated file: $KEEPSAMPLES (POP=$POP)"
echo "[INFO] Base PFILE: $PFILE"
echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] SNPs-only: $SNPS_ONLY"
echo "[INFO] Drop ambiguous: $DROP_AMBIG"
echo "[INFO] HWE command: $HWE_COMMAND"
echo "[INFO] Geno threshold: $GENO_THRES"
echo "[INFO] MAF threshold: $MAF_THRES"
echo "[INFO] Threads: $THREADS"
echo "[INFO] VCF CHR prefix: $VCF_CHR_PREFIX"

# Step 1: Subset samples
./bin/plink2 --pfile "$PFILE" --keep "$KEEPSAMPLES" \
  --make-pgen --threads "$THREADS" --out "$OUTDIR/kept"

INPFX="$OUTDIR/kept"

# Step 2: SNPs-only filter
if [[ "$SNPS_ONLY" == "yes" ]]; then
  ./bin/plink2 --pfile "$INPFX" --snps-only just-acgt --max-alleles 2 \
    --make-pgen --threads "$THREADS" --out "${INPFX}_snps"
  INPFX="${INPFX}_snps"
fi

# Step 3: Drop ambiguous SNPs
if [[ "$DROP_AMBIG" == "yes" ]]; then
  awk 'BEGIN{FS=OFS="\t"}
    $1 !~ /^#/ {
      ref=toupper($4); alt=toupper($5)
      if (alt ~ /,/) next
      if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) next
      amb = ((ref=="A"&&alt=="T")||(ref=="T"&&alt=="A")||(ref=="C"&&alt=="G")||(ref=="G"&&alt=="C"))
      if (!amb) print $3
    }' "$INPFX.pvar" > "${INPFX}_nonambig.ids"
  ./bin/plink2 --pfile "$INPFX" --extract "${INPFX}_nonambig.ids" \
    --make-pgen --threads "$THREADS" --out "${INPFX}_noambig"
  INPFX="${INPFX}_noambig"
fi

# Step 4: HWE
if [[ -n "$HWE_COMMAND" ]]; then
  ./bin/plink2 --pfile "$INPFX" --hwe $HWE_COMMAND --mac 20 --threads "$THREADS" \
    --write-snplist --out "${INPFX}_keep_hwe"

  ./bin/plink2 --pfile "$INPFX" --extract "${INPFX}_keep_hwe.snplist" \
    --make-pgen --threads "$THREADS" --out "${INPFX}_hwe"
  INPFX="${INPFX}_hwe"
fi

# Step 5: Variant QC filters + sample QC with good variants
PLINK_ARGS=(--pfile "$INPFX" --make-pgen --mind 0.02 --sort-vars --threads "$THREADS" --out "${INPFX}_qc")
[[ -n "$FA" ]] && PLINK_ARGS+=(--ref-from-fa force --fa "$FA")
[[ -n "$GENO_THRES" ]] && PLINK_ARGS+=(--geno "$GENO_THRES")
[[ -n "$MAF_THRES" ]] && PLINK_ARGS+=(--maf "$MAF_THRES")

./bin/plink2 "${PLINK_ARGS[@]}"

echo "[INFO] Finished $POP QC: ${INPFX}_qc.{pgen,pvar,psam}"

if [[ "$VCF_OUT" == "yes" ]]; then
  # rename the variant IDs to a standard format for imputation servers
  ./bin/plink2 --pfile "${INPFX}_qc" \
    --set-all-var-ids 'chr@:#:$r:$a' \
    --make-pgen --threads "$THREADS" --out "${INPFX}_qc_renamed"
  echo "[INFO] Exporting final QCed data to VCF"
  mkdir -p "${OUTDIR}/vcf"
  # Get available chromosomes from .pvar file (skip header lines, handle chr prefix)
  CHRS=$(awk 'NR>1 && $1 !~ /^#/ {c=$1; sub(/^chr/,"",c); print c}' "${INPFX}_qc.pvar" | sort -u)
  echo "[INFO] Exporting chromosomes: $CHRS"
  for c in $CHRS; do
    ./bin/plink2 --pfile "${INPFX}_qc_renamed" --chr "$c" --output-chr chrM --recode vcf bgz \
      --out "${OUTDIR}/vcf/chr${c}"
  done
  awk 'NR>1{print $1,$2}' OFS='\t' "${INPFX}_qc_renamed.psam" > "${OUTDIR}/vcf/sample_info.list"
fi


echo -e "\n==> DONE. QC outputs in $OUTDIR/*_qc/"
