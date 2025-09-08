#!/usr/bin/env bash
set -euo pipefail

# Merge hg38-normalized study data with a reference panel and compute PCs.
# Requirements: ./bin/plink2 (new build), both datasets in hg38 and normalized;
#               variant IDs as chr:pos:ref:alt (or you let this script set them).

# Inputs
# --study-pfile   PATH   (required) study dataset prefix (PGEN/PVAR/PSAM) — hg38
# --ref-bfile     PATH   (required) reference dataset prefix (BED/BIM/FAM) — already hg38; IDs ideally chr:pos:ref:alt
# --fa            PATH   (optional) hg38 FASTA (bgzip+.fai) to re-normalize/align both sides defensively
# --drop-ambig    yes|no (default yes) drop ambiguous A/T & C/G SNPs before merging
# --threads       INT    (default 2)
# --outdir        PATH   (default merge_ref_out)
# --out-prefix    STR    (default merged_study_ref)
#
# Output: ${outdir}/${outprefix}.pgen/.pvar/.psam (merged),
#         plus LD-pruned set and PCA outputs.

STUDY=""; REF=""; FA=""
DROP_AMBIG="yes"
THREADS=2
OUTDIR="merge_ref_out"
OUTPFX="merged_study_ref"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --study-pfile) STUDY="$2"; shift 2 ;;
    --ref-bfile)   REF="$2";   shift 2 ;;
    --fa)          FA="$2";    shift 2 ;;
    --drop-ambig)  DROP_AMBIG="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --outdir)      OUTDIR="$2"; shift 2 ;;
    --out-prefix)  OUTPFX="$2"; shift 2 ;;
    -h|--help) sed -n '1,200p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -n "$STUDY" ]] || { echo "ERROR: --study-pfile required"; exit 1; }
[[ -n "$REF"   ]] || { echo "ERROR: --ref-bfile required"; exit 1; }
if [[ -n "$FA" ]]; then
  [[ -f "$FA" ]] || { echo "ERROR: FASTA not found: $FA"; exit 1; }
  [[ -f "$FA.fai" ]] || { echo "ERROR: FASTA index missing: $FA.fai"; exit 1; }
fi
command -v ./bin/plink2 >/dev/null || { echo "ERROR: ./bin/plink2 not in PATH"; exit 1; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/merge.log"; exec > >(tee -a "$LOG") 2>&1
trap 'echo "[ERROR] Failed at line $LINENO" >&2' ERR

echo "[INFO] Merge start"
echo "[INFO] study=$STUDY  ref(bfile)=$REF  threads=$THREADS  drop_ambig=$DROP_AMBIG"

# 0) Convert reference bfile -> pgen
./bin/plink2 --bfile "$REF" --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_pgen"

# 1) Defensively normalize/align + standardize IDs on both sides
#    (Safe even if already normalized; ensures REF/ALT + IDs consistent)
if [[ -n "$FA" ]]; then
  ./bin/plink2 --pfile "$STUDY" \
         --fa "$FA" --normalize --ref-from-fa force \
         --autosome \
         --snps-only just-acgt \
         --maf 0.1 --geno 0.01 \
         --set-all-var-ids 'chr@:#:$r:$a' \
         --sort-vars --make-pgen --threads "$THREADS" --out "$OUTDIR/study_norm"

  ./bin/plink2 --pfile "$OUTDIR/ref_pgen" \
         --fa "$FA" --normalize --ref-from-fa force \
         --autosome \
         --snps-only just-acgt \
         --maf 0.1 --geno 0.01 \
         --set-all-var-ids 'chr@:#:$r:$a' \
         --sort-vars --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_norm"
  STUDY="$OUTDIR/study_norm"
  REF="$OUTDIR/ref_norm"
else
  # At least standardize IDs for both
  ./bin/plink2 --pfile "$OUTDIR/ref_pgen" \
         --set-all-var-ids 'chr@:#:$r:$a' --sort-vars \
         --autosome --snps-only just-acgt \
         --maf 0.1 --geno 0.01 \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_named"
  ./bin/plink2 --pfile "$STUDY" \
         --set-all-var-ids 'chr@:#:$r:$a' --sort-vars \
         --autosome --snps-only just-acgt \
         --maf 0.1 --geno 0.1 \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/study_named"
  STUDY="$OUTDIR/study_named"
  REF="$OUTDIR/ref_named"
fi

# 2) Optional: drop ambiguous A/T & C/G SNPs (helps if any strand quirks remain)
if [[ "$DROP_AMBIG" == "yes" ]]; then
  awk 'BEGIN{FS=OFS="\t"}
     $1 !~ /^#/ {
       ref=toupper($4); alt=toupper($5)
       # skip multi-allelic and non-ACGT
       if (alt ~ /,/) next
       if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) next
       amb = ((ref=="A"&&alt=="T")||(ref=="T"&&alt=="A")||(ref=="C"&&alt=="G")||(ref=="G"&&alt=="C"))
       if (!amb) print $3
     }' "$STUDY.pvar" > "$OUTDIR/keep_nonambig.ids"
  ./bin/plink2 --pfile "$STUDY" --extract "$OUTDIR/keep_nonambig.ids" \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/study_noambig"
  ./bin/plink2 --pfile "$REF" --extract "$OUTDIR/keep_nonambig.ids" \
         --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_noambig"
  STUDY="$OUTDIR/study_noambig"
  REF="$OUTDIR/ref_noambig"
fi

# 3) Subset study by ref.ids
cut -f3 "$REF.pvar"   | sed '1d' | sort > "$OUTDIR/ref.ids"
./bin/plink2 --pfile "$STUDY"\
       --extract "$OUTDIR/ref.ids" \
       --sort-vars --make-pgen --threads "$THREADS" --out "$OUTDIR/study_shared"
./bin/plink2 --pfile "$OUTDIR/study_shared" --make-bed --threads "$THREADS" --out "$OUTDIR/study_shared_b"  # also make bed for merging step
# 4) Subset ref by shared.ids
grep -v '^#' "$OUTDIR/study_shared.pvar" | cut -f3 > "$OUTDIR/shared.ids"
./bin/plink2 --pfile "$REF" \
       --extract "$OUTDIR/shared.ids" \
       --sort-vars --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_shared"
./bin/plink2 --pfile "$OUTDIR/ref_shared" --make-bed --threads "$THREADS" --out "$OUTDIR/ref_shared_b"  # also make bed for merging step

NSHARED=$(wc -l < "$OUTDIR/shared.ids")
echo "[INFO] Shared variants: $NSHARED"

# 5) Merge (stack samples) with --pmerge-list
#    Base = reference; pmerge the study on top (or vice versa—either is fine)
echo "$OUTDIR/study_shared" > "$OUTDIR/pmerge_list.txt"
plink --bfile "$OUTDIR/ref_shared_b" \
       --bmerge "$OUTDIR/study_shared_b" \
       --make-bed --threads "$THREADS" --out "$OUTDIR/$OUTPFX"

# 7) LD prune on merged data (autosomes + X is typical; adjust if needed)
./bin/plink2 --bfile "$OUTDIR/$OUTPFX" \
       --autosome \
       --maf 0.1 --geno 0.01 \
       --indep-pairwise 200 50 0.2 \
       --threads "$THREADS" --out "$OUTDIR/$OUTPFX.prune"

# 8) PCA on merged, pruned set
N_MERGED=$(wc -l < "$OUTDIR/$OUTPFX.fam")
PCA_FLAG=()
if [[ "$N_MERGED" -gt 5000 ]]; then PCA_FLAG=(approx); echo "[INFO] PCA: using approx (N=$N_MERGED > 5000)"; else echo "[INFO] PCA: using exact (N=$N_MERGED ≤ 5000)"; fi

./bin/plink2 --bfile "$OUTDIR/$OUTPFX" \
       --extract "$OUTDIR/$OUTPFX.prune.prune.in" \
       --pca 20 ${PCA_FLAG:+${PCA_FLAG[@]}} \
       --threads "$THREADS" --out "$OUTDIR/$OUTPFX.pca"


echo -e "\n==> DONE. Outputs in $OUTDIR\n- Merged PGEN: $OUTDIR/$OUTPFX.{pgen,pvar,psam}\n- Prune list:  $OUTDIR/$OUTPFX.prune.prune.in\n- PCA:         $OUTDIR/$OUTPFX.pca.{eigenvec,eigenval}"
