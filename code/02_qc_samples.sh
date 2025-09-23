# 02_qc_samples.sh — ./bin/plink2 sample QC → final NON-PRUNED dataset (no PCA)
#!/usr/bin/env bash
set -euo pipefail

# Inputs:
# --pfile PATH          (required) prepped hg38 dataset from 01_prepare_variants.sh
# --fa PATH             (optional) hg38 FASTA (bgzip+.fai); used only when sex-check runs
# --keep-samples PATH   (optional) sample list to retain (IID-only or "FID IID")
# --mind FLOAT          (default 0.02) sample missingness
# --indep-kb INT        (default 200)
# --indep-step INT      (default 50)
# --indep-r2 FLOAT      (default 0.2)
# --het-z FLOAT         (default 4.0) |Z| cutoff for F coefficient outliers
# --king-cutoff FLOAT   (default 0.0884) (~3rd-degree) for unrelated set
# --threads N           (default 2)
# --outdir DIR          (default qc_samples_out))

PFILE=""; FA=""; SKEEP=""
MIND=0.02
INDEP_KB=200; INDEP_STEP=50; INDEP_R2=0.2
HET_Z=4.0
KING=0.0884
THREADS=2
OUTDIR="qc_samples_out"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfile) PFILE="$2"; shift 2 ;;
    --fa) FA="$2"; shift 2 ;;
    --keep-samples) SKEEP="$2"; shift 2 ;;
    --mind) MIND="$2"; shift 2 ;;
    --indep-kb) INDEP_KB="$2"; shift 2 ;;
    --indep-step) INDEP_STEP="$2"; shift 2 ;;
    --indep-r2) INDEP_R2="$2"; shift 2 ;;
    --het-z) HET_Z="$2"; shift 2 ;;
    --king-cutoff) KING="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) sed -n '1,200p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$PFILE" ]] && { echo "ERROR: --pfile required"; exit 1; }
if [[ -n "$FA" ]]; then
  [[ -f "$FA" ]] || { echo "ERROR: FASTA not found: $FA"; exit 1; }
  [[ -f "$FA.fai" ]] || { echo "ERROR: FASTA index missing: $FA.fai"; exit 1; }
fi
if [[ -n "$SKEEP" ]]; then
  [[ -f "$SKEEP" ]] || { echo "ERROR: keep-samples file not found: $SKEEP"; exit 1; }
fi
command -v ./bin/plink2 >/dev/null || { echo "ERROR: ./bin/plink2 not in PATH (expected local binary)"; exit 1; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/qc_samples.log"; exec > >(tee -a "$LOG") 2>&1
trap 'echo "[ERROR] Failed at line $LINENO" >&2' ERR

echo "[INFO] Sample QC start"
echo "[INFO] Input=$PFILE  Threads=$THREADS  Outdir=$OUTDIR  mind=$MIND  het_z=$HET_Z  king=$KING"

# Clean any old lists
rm -f "$OUTDIR"/_remove_* "$OUTDIR"/final_remove.list "$OUTDIR"/final_keep.list

# --------------------------------------------------------------------
# 0) Optional initial keep filter (IID-only or FID IID)
INPFX="$PFILE"
if [[ -n "$SKEEP" ]]; then
  ./bin/plink2 --pfile "$INPFX" --keep "$SKEEP" \
               --make-pgen --threads "$THREADS" --out "$OUTDIR/00_kept"
  INPFX="$OUTDIR/00_kept"
fi
BASE_FULL="$INPFX"   # BASE_FULL = starting point for all sample-QC steps

# 1) Sample call rate (mind) — this defines the FULL baseline for later removals
./bin/plink2 --pfile "$INPFX" --mind "$MIND" \
             --make-pgen --threads "$THREADS" --out "$OUTDIR/01_mind"
INPFX="$OUTDIR/01_mind"
# Build mind removal list
awk -v label="mind_$MIND" 'NR==1{next} {print $1, $2, label}' OFS='\t' "$INPFX.mindrem.id" \
  | sort -k1,1n -k2,2n > "$OUTDIR/_remove_mind.list" || true




# --------------------------------------------------------------------
# 2) Sex check (only if chrX present and SEX usable) → produce removal list
HAS_CHRX=0
if [[ -f "$INPFX.pvar" ]]; then
  if awk 'BEGIN{FS="\t"} $1!~/^#/ { if($1=="chrX"){found=1; exit} } END{exit !found}' "$INPFX.pvar"; then
    HAS_CHRX=1
  fi
fi

HAS_SEX=0
if [[ -f "$INPFX.psam" ]]; then
  SEXCOL=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SEX"){print i; exit}}' "$INPFX.psam" || true)
  if [[ -n "$SEXCOL" ]]; then
    NSEX=$(awk -v c="$SEXCOL" 'NR>1 && $c!=0{n++} END{print n+0}' "$INPFX.psam" 2>/dev/null || echo 0)
    [[ "$NSEX" -gt 0 ]] && HAS_SEX=1
  fi
fi

if [[ "$HAS_CHRX" -eq 1 && "$HAS_SEX" -eq 1 ]]; then
  echo "[INFO] chrX & SEX present → chrX-only with PAR split, then sex check"
  ./bin/plink2 --pfile "$INPFX" --split-par hg38 --chr X \
               --make-pgen --threads "$THREADS" --out "$OUTDIR/02_chrx_split"

  if [[ -n "${FA:-}" ]]; then
    ./bin/plink2 --pfile "$OUTDIR/02_chrx_split" --check-sex --fa "$FA" \
                 --threads "$THREADS" --out "$OUTDIR/02_sexcheck"
  else
    ./bin/plink2 --pfile "$OUTDIR/02_chrx_split" --check-sex \
                 --threads "$THREADS" --out "$OUTDIR/02_sexcheck"
  fi

  if [[ -s "$OUTDIR/02_sexcheck.sexcheck" ]]; then
    awk 'NR>1 && $NF=="PROBLEM"{print $1}' "$OUTDIR/02_sexcheck.sexcheck" \
      | sort -u > "$OUTDIR/_remove_sex_fail.list" || true
  fi
else
  [[ "$HAS_CHRX" -ne 1 ]] && echo "[INFO] chrX not present → skipping sex check"
  [[ "$HAS_SEX"  -ne 1 ]] && echo "[INFO] SEX not present/usable → skipping sex check"
fi

# --------------------------------------------------------------------
# 3) LD pruning (for het/rel only) — pruned set is for stats, not for final data
./bin/plink2 --pfile "$INPFX" \
             --autosome \
             --geno 0.01 \
             --maf 0.05 \
             --indep-pairwise "$INDEP_KB" "$INDEP_STEP" "$INDEP_R2" \
             --threads "$THREADS" --out "$OUTDIR/03_prune"

./bin/plink2 --pfile "$INPFX" --extract "$OUTDIR/03_prune.prune.in" \
             --make-pgen --threads "$THREADS" --out "$OUTDIR/03_pruned"

# --------------------------------------------------------------------
# 4) Heterozygosity on pruned set → flag outliers (build list only)
./bin/plink2 --pfile "$OUTDIR/03_pruned" --het \
             --threads "$THREADS" --out "$OUTDIR/04_het"

awk -v z="$HET_Z" '
  BEGIN{ FS="[ \t]+"; OFS="\t" }
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="#FID") fidc=i;
      if($i=="IID") iidc=i;
      if($i=="F")   fc=i;
    }
    next
  }
  {
    f = $(fc)
    if(f!="NA" && f!="nan"){
      n++; sum+=f; s2+=f*f
      key = $(fidc) SUBSEP $(iidc)
      vals[key] = f
    }
  }
  END{
    if(n>1){
      mean = sum/n
      sd = sqrt((s2 - n*mean*mean)/(n-1))
      for(k in vals){
        f = vals[k]
        zval = (sd>0)? (f-mean)/sd : 0
        if(zval>=z || zval<=-z){
          split(k, a, SUBSEP)   # a[1]=FID, a[2]=IID
          print a[1], a[2], "het_Z_" z
        }
      }
    }
  }
' "$OUTDIR/04_het.het" | LC_ALL=C sort -t $'\t' -k1,1n -k2,2n > "$OUTDIR/_remove_het_outlier.list" || true

# --------------------------------------------------------------------
# 5) Relatedness (KING) on pruned set → keep unrelated; relateds = current minus unrelated
./bin/plink2 --pfile "$OUTDIR/03_pruned" --king-cutoff "$KING" \
             --make-pgen --threads "$THREADS" --out "$OUTDIR/05_unrelated"

# Build related-remove list
awk -v label="king_$KING" 'NR==1{next} {print $1, $2, label}' OFS='\t' "$OUTDIR/05_unrelated.king.cutoff.out.id" \
  | sort -k1,1n -k2,2n > "$OUTDIR/_remove_related.list" || true

# --------------------------------------------------------------------
# 6) Combine remove lists → final_remove.list (unique, sorted)
cat \
  "$OUTDIR"/_remove_mind.list \
  "$OUTDIR"/_remove_sex_fail.list \
  "$OUTDIR"/_remove_het_outlier.list \
  "$OUTDIR"/_remove_related.list \
  2>/dev/null > "$OUTDIR/final_remove.list" || true

N_MIND=$(wc -l < "$OUTDIR/_remove_mind.list" 2>/dev/null || echo 0)
N_SEX=$(wc -l < "$OUTDIR/_remove_sex_fail.list" 2>/dev/null || echo 0)
N_HET=$(wc -l < "$OUTDIR/_remove_het_outlier.list" 2>/dev/null || echo 0)
N_REL=$(wc -l < "$OUTDIR/_remove_related.list" 2>/dev/null || echo 0)
N_REM=$(wc -l < "$OUTDIR/final_remove.list" 2>/dev/null || echo 0)

awk -F'\t' '{
  if (NF == 3) {
    # FID IID reason
    print $1, $2
  } else if (NF == 2) {
    # IID reason
    print $1
  }
}' "$OUTDIR/final_remove.list" > "$OUTDIR/final_remove.id"


# --------------------------------------------------------------------
# 7) Produce FINAL NON-PRUNED QC dataset by removing final_remove from BASE_FULL
if [[ -s "$OUTDIR/final_remove.id" ]]; then
  ./bin/plink2 --pfile "$PFILE" --remove "$OUTDIR/final_remove.list" \
               --make-pgen --threads "$THREADS" --out "$OUTDIR/final_keep"
else
  # No removals → just copy BASE_FULL forward
  ./bin/plink2 --pfile "$PFILE" \
               --make-pgen --threads "$THREADS" --out "$OUTDIR/final_keep"
fi

# --------------------------------------------------------------------
# Summary
NS_START=$(($(wc -l < "$BASE_FULL.psam")-1))
NS_FINAL=$(($(wc -l < "$OUTDIR/final_keep.psam")-1))

echo "[INFO] N samples at start (after keep if any): $NS_START"
echo "[INFO] Removals: mind=$N_MIND"
echo "[INFO] Removals: sex=$N_SEX"
echo "[INFO] Removals: het=$N_HET"
echo "[INFO] Removals: related=$N_REL"
echo "[INFO] Removals: total=$N_REM"
echo "[INFO] N samples after removal: $NS_FINAL"
echo "[INFO] Final PGEN: $OUTDIR/final_keep.{pgen,pvar,psam}"
echo "[INFO] Remove list: $OUTDIR/final_remove.list"
