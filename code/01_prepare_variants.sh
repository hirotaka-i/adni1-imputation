# 01_prepare_variants.sh — ./bin/plink2-first sanity + (optional) liftOver + normalize (hg38)
#!/usr/bin/env bash
set -euo pipefail

# Inputs:
# --bfile / --pfile     (one required)
# --fa PATH             (bgzip + .fai; hg38 target FASTA)
# --chain PATH          (optional; e.g. hg19ToHg38.over.chain.gz)
# --chain-back PATH     (optional; e.g. hg38ToHg19.over.chain.gz; if doing two-way liftover)
# --vartokeep PATH      (optional; variant IDs to keep; one per line)
# --keep-samples PATH   (optional; sample list to keep: IID only (preferred) or FID IID)
# --geno-pre FLOAT      (default 0.10) early per-variant missingness
# --dup-policy STR      (force-first|exclude-all) default force-first
# --threads N           (default 2)
# --outdir DIR          (default prep_vars_out)

BFILE=""; PFILE=""; FA=""; CHAIN=""; CHAIN_BACK=""; VKEEP=""; SKEEP=""
THREADS=2; OUTDIR="prep_vars_out"; GENO_PRE=0.10; DUP_POLICY="force-first"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) BFILE="$2"; shift 2 ;;
    --pfile) PFILE="$2"; shift 2 ;;
    --fa) FA="$2"; shift 2 ;;
    --chain) CHAIN="$2"; shift 2 ;;
    --chain-back) CHAIN_BACK="$2"; shift 2 ;;
    --vartokeep) VKEEP="$2"; shift 2 ;;
    --keep-samples) SKEEP="$2"; shift 2 ;;
    --geno-pre) GENO_PRE="$2"; shift 2 ;;
    --dup-policy) DUP_POLICY="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) sed -n '1,200p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

# --- sanity checks ---
[[ -z "$FA" ]] && { echo "ERROR: --fa required (bgzip + .fai)"; exit 1; }
[[ -z "$BFILE" && -z "$PFILE" ]] && { echo "ERROR: provide --bfile or --pfile"; exit 1; }
[[ -f "$FA" ]] || { echo "ERROR: FASTA not found: $FA"; exit 1; }
[[ -f "$FA.fai" ]] || { echo "ERROR: FASTA index missing: $FA.fai"; exit 1; }
if [[ -n "$CHAIN" ]]; then [[ -f "$CHAIN" ]] || { echo "ERROR: chain not found: $CHAIN"; exit 1; }; fi
if [[ -n "$CHAIN_BACK" ]]; then [[ -f "$CHAIN_BACK" ]] || { echo "ERROR: chain-back not found: $CHAIN_BACK"; exit 1; }; fi
if [[ -n "$VKEEP" ]]; then [[ -f "$VKEEP" ]] || { echo "ERROR: --vartokeep not found: $VKEEP"; exit 1; }; fi
if [[ -n "$SKEEP" ]]; then [[ -f "$SKEEP" ]] || { echo "ERROR: --keep-samples not found: $SKEEP"; exit 1; }; fi
case "$DUP_POLICY" in force-first|exclude-all) ;; *) echo "ERROR: --dup-policy must be force-first|exclude-all"; exit 1;; esac

mkdir -p "$OUTDIR"
LOG="$OUTDIR/prep_variants.log"; exec > >(tee -a "$LOG") 2>&1
trap 'echo "[ERROR] Failed at line $LINENO" >&2' ERR

inflag=(); [[ -n "$BFILE" ]] && inflag+=(--bfile "$BFILE") || inflag+=(--pfile "$PFILE")

echo "[INFO] Threads=$THREADS  Outdir=$OUTDIR  Geno-pre=$GENO_PRE  Dup-policy=$DUP_POLICY"
[[ -n "$CHAIN" ]] && echo "[INFO] Liftover chain: $CHAIN" || echo "[INFO] Liftover: SKIPPED"
[[ -n "$SKEEP" ]] && echo "[INFO] Keeping samples listed in: $SKEEP"

# 0) Initial reports (pre-filters, whole dataset)
./bin/plink2 "${inflag[@]}" --freq --missing --hardy --threads "$THREADS" --out "$OUTDIR/00_initial"

# 1) Variant de-dup (policy) + write PGEN
./bin/plink2 "${inflag[@]}" \
  --rm-dup "$DUP_POLICY" \
  --make-pgen --threads "$THREADS" --out "$OUTDIR/01_nodup"

# 1b) Optional sample keep (IID-only supported by ./bin/plink2; FID IID also OK)
INPFX="$OUTDIR/01_nodup"
if [[ -n "$SKEEP" ]]; then
  ./bin/plink2 --pfile "$INPFX" \
         --keep "$SKEEP" \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/01b_keep"
  INPFX="$OUTDIR/01b_keep"
fi

# 2) Optional variant include list
if [[ -n "$VKEEP" ]]; then
  ./bin/plink2 --pfile "$INPFX" --extract "$VKEEP" \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/01c_varfilt"
  INPFX="$OUTDIR/01c_varfilt"
fi

# 3) Lenient per-variant missingness (only biallelic snps retained)
./bin/plink2 --pfile "$INPFX" --geno "$GENO_PRE" \
       --snps-only just-acgt --max-alleles 2 \
       --make-pgen --threads "$THREADS" --out "$OUTDIR/02_lenient"

# 4) Optional liftover via UCSC liftOver (BED: 0-based; use full REF span)
LIFTEDPFX="$OUTDIR/02_lenient"
if [[ -n "$CHAIN" ]]; then
  echo "[INFO] Running liftOver"
  ./bin/plink2 --pfile "$LIFTEDPFX" \
         --set-all-var-ids '@:#:$r:$a' \
         --new-id-max-allele-len 999 truncate \
         --rm-dup "$DUP_POLICY" \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/03_tmp_named"

  # BED input for liftOver: add 'chr' for UCSC chains, length-aware for indels
  awk 'NR>1{
    chr=$1; pos=$2; id=$3; ref=$4;
    if (chr !~ /^chr/) chr="chr" chr;
    start=pos-1; end=start+length(ref);
    print chr, start, end, id
  }' OFS='\t' "$OUTDIR/03_tmp_named.pvar" > "$OUTDIR/03_lift.in.bed"

  command -v liftOver >/dev/null || { echo "ERROR: liftOver not in PATH"; exit 1; }
  liftOver "$OUTDIR/03_lift.in.bed" "$CHAIN" \
           "$OUTDIR/03_lift.out.bed" "$OUTDIR/03_lift.unmapped"

  # Dedup potential multi-maps per ID (keep first)
  awk '!seen[$4]++' "$OUTDIR/03_lift.out.bed" > "$OUTDIR/03_lift.out.unique.bed"

  echo "[INFO] liftover counts: in=$(wc -l < "$OUTDIR/03_lift.in.bed")  out=$(wc -l < "$OUTDIR/03_lift.out.unique.bed")  unmapped=$(wc -l < "$OUTDIR/03_lift.unmapped")"
  [[ -s "$OUTDIR/03_lift.out.unique.bed" ]] || { echo "ERROR: liftover produced 0 mappings (check chain/chr style)"; exit 1; }
  BED="$OUTDIR/03_lift.out.unique.bed"

  # chain back (optional)
  if [[ -n "$CHAIN_BACK" ]]; then
       echo "[INFO] Running liftOver back-mapping"
       liftOver "$OUTDIR/03_lift.out.unique.bed" "$CHAIN_BACK" "$OUTDIR/03_lift.back.bed" "$OUTDIR/03_lift.back.unmapped"
       # build ID->cood maps for original (pre-liftover input and back-lifted)
       awk '{print $4, $1, $2, $3}' OFS='\t' "$OUTDIR/03_lift.in.bed" | sort -k1,1 > "$OUTDIR/03_lift.orig.map"
       awk '{print $4, $1, $2, $3}' OFS='\t' "$OUTDIR/03_lift.back.bed" | sort -k1,1 > "$OUTDIR/03_lift.back.map"
       # compare original vs back-lifted positions; keep only perfect matches
       join -t $'\t' -1 1 -2 1 "$OUTDIR/03_lift.orig.map" "$OUTDIR/03_lift.back.map" \
              | awk -F'\t' '$2==$5 && $3==$6 && $4==$7 {print $1}' > "$OUTDIR/03_lift.roundtrip.ids"
       echo "[INFO] Liftover back-mapping matches: $(wc -l < "$OUTDIR/03_lift.roundtrip.ids")"
       # filter to only those that perfectly back-mapped
       awk 'NR==FNR{keep[$1]=1; next} $4 in keep' "$OUTDIR/03_lift.roundtrip.ids" "$OUTDIR/03_lift.out.unique.bed" \
              > "$OUTDIR/03_lift.out.unique.filtered.bed"
       echo "[INFO] After back-mapping filter: $(wc -l < "$OUTDIR/03_lift.out.unique.filtered.bed") variants remain"
       BED="$OUTDIR/03_lift.out.unique.filtered.bed"
  fi


  cut -f4 "$BED" > "$OUTDIR/03_keep.ids"
  awk '{print $4, $2+1}' OFS='\t' "$BED" > "$OUTDIR/03_update.pos"
  # Create chr update map: variantID -> newChr (keep chr prefix from liftOver output)
  awk '{print $4, $1}' OFS='\t' "$BED" > "$OUTDIR/03_update.chr"

  ./bin/plink2 --pfile "$OUTDIR/03_tmp_named" \
         --extract "$OUTDIR/03_keep.ids" \
         --update-map "$OUTDIR/03_update.pos" \
         --update-chr "$OUTDIR/03_update.chr" 2 1 \
         --sort-vars \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/04_hg38"

  LIFTEDPFX="$OUTDIR/04_hg38"
else
  # No liftover: add chr prefix directly since data likely doesn't have it
  echo "[INFO] No liftover - adding 'chr' prefix to match FASTA reference"
  awk 'NR==1 {print; next} {print $1, "chr"$1}' "$LIFTEDPFX.pvar" > "$OUTDIR/04_chr_update.txt"
  ./bin/plink2 --pfile "$LIFTEDPFX" \
         --update-chr "$OUTDIR/04_chr_update.txt" \
         --sort-vars \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/04_hg38"
  LIFTEDPFX="$OUTDIR/04_hg38"
fi

# At this point, LIFTEDPFX should have chr-prefixed chromosomes matching the FASTA

# 5) Normalize vs hg38 FASTA, align REF/ALT, dedup (policy)
# NOTE: ref-from-fa only aligns to FASTA, NOT to imputation panel strand
# Use HRC-1000G checking (script 05) after this for proper panel alignment
echo "[INFO] Normalizing against FASTA reference (chr-prefixed chromosomes)"
./bin/plink2 --pfile "$LIFTEDPFX" --fa "$FA" \
       --normalize --ref-from-fa force \
       --rm-dup "$DUP_POLICY" --sort-vars \
       --make-pgen --threads "$THREADS" --out "$OUTDIR/05_hg38_norm_aligned"

# 6) Set variant IDs
./bin/plink2 --pfile "$OUTDIR/05_hg38_norm_aligned" \
       --set-all-var-ids 'chr@:#:$r:$a' \
       --new-id-max-allele-len 999 truncate \
       --make-pgen --threads "$THREADS" --out "$OUTDIR/06_hg38_named"

# 6) Final dup sweep with chosen policy
./bin/plink2 --pfile "$OUTDIR/06_hg38_named" \
       --rm-dup "$DUP_POLICY" \
       --chr 1-22,X,Y,MT \
       --make-pgen --threads "$THREADS" --out "$OUTDIR/hg38_prepped"

echo -e "\n==> DONE.\nReport: $LOG\nClean PGEN: $OUTDIR/hg38_prepped.{pgen,pvar,psam}\n"
