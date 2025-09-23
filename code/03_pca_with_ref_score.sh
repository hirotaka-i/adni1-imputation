#!/usr/bin/env bash
set -euo pipefail

# Train PCA on a reference panel (REF) and project study samples (STUDY),
# producing PCs on the SAME numeric scale for both cohorts.
#
# Pipeline summary (unchanged from your working version):
#   0) Convert REF bfile -> PGEN for harmonization steps
#   1) Normalize/align both to FASTA (if provided), autosomes, A/C/G/T, MAF/GENO, set IDs, sort
#   2) (Optional) Drop ambiguous palindromic SNPs
#   3) Intersect IDs, align alleles, create shared BEDs
#   3.5) Merge to obtain a "good by merge" SNP list
#   4) LD-prune on REF only; apply same prune list to STUDY
#   5) PLINK2 PCA on REF with var-wts; project REF and STUDY via --score variance-standardize
#   6) Assemble combined outputs for plotting
#
# Requirements:
#   - ./bin/plink2 (v2) and plink (v1.9+)
#   - REF and STUDY are hg38; pass --fa to re-normalize/force REF/ALT if desired
# Inputs
#   --study-pfile   PATH   (required) study dataset prefix (PGEN/PVAR/PSAM) — hg38
#   --ref-bfile     PATH   (required) reference dataset prefix (BED/BIM/FAM) — hg38
#   --fa            PATH   (optional) hg38 FASTA (bgzip+.f
#   --drop-ambig   STR    (default yes) drop ambiguous A/T & C/G SNPs? yes|no
#   --threads      N      (default 2)
#   --outdir       DIR    (default pca_proj_out)
#   --out-prefix   STR    (default study_vs_ref)
#   --k            INT    (default 20) number of PCs to compute
#   --prune-win    INT    (default 200) LD prune window size (
#   --prune-step   INT    (default 50) LD prune step size
#   --prune-r2     FLOAT  (default 0.2) LD prune
#   --exclude-range PATH (optional) exclude regions (e.g. inversions) from LD prune
#   --maf          FLOAT  (default 0.10) MAF filter (both
#   --geno         FLOAT  (default 0.01) per-variant missingness (both)

# --------- CLI ---------
STUDY=""; REF=""; FA=""
DROP_AMBIG="yes"
THREADS=2
OUTDIR="pca_proj_out"
OUTPFX="study_vs_ref"
K=20
PRUNE_WIN=200
PRUNE_STEP=50
PRUNE_R2=0.2
EXCL_RANGE=""
MAF=0.10
GENO=0.01

while [[ $# -gt 0 ]]; do
  case "$1" in
    --study-pfile) STUDY="$2"; shift 2 ;;
    --ref-bfile)   REF="$2";   shift 2 ;;
    --fa)          FA="$2";    shift 2 ;;
    --drop-ambig)  DROP_AMBIG="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --outdir)      OUTDIR="$2"; shift 2 ;;
    --out-prefix)  OUTPFX="$2"; shift 2 ;;
    --k)           K="$2"; shift 2 ;;
    --prune-win)   PRUNE_WIN="$2"; shift 2 ;;
    --prune-step)  PRUNE_STEP="$2"; shift 2 ;;
    --prune-r2)    PRUNE_R2="$2"; shift 2 ;;
    --exclude-range) EXCL_RANGE="$2"; shift 2 ;;
    --maf)         MAF="$2"; shift 2 ;;
    --geno)        GENO="$2"; shift 2 ;;
    -h|--help) sed -n '1,200p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -n "$STUDY" ]] || { echo "ERROR: --study-pfile required"; exit 1; }
[[ -n "$REF"   ]] || { echo "ERROR: --ref-bfile required"; exit 1; }
if [[ -n "$FA" ]]; then
  [[ -f "$FA" && -f "$FA.fai" ]] || { echo "ERROR: FASTA/FAI missing"; exit 1; }
fi
command -v ./bin/plink2 >/dev/null || { echo "ERROR: ./bin/plink2 not found"; exit 1; }
command -v plink >/dev/null || { echo "ERROR: plink (1.9+) not found"; exit 1; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/flashpca_train_project.log"; exec > >(tee -a "$LOG") 2>&1
trap 'echo "[ERROR] Failed at line $LINENO" >&2' ERR
info(){ echo "[INFO] $*"; }

info "Start PCA train→project"
info "study=$STUDY  ref(bfile)=$REF  threads=$THREADS  K=$K  maf=$MAF  geno=$GENO  drop_ambig=$DROP_AMBIG"

# 0) REF bfile -> PGEN
./bin/plink2 --bfile "$REF" --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_pgen"

# 1) Normalize/QC/IDs/sort (symmetric)
if [[ -n "$FA" ]]; then
  info "Normalizing with FASTA; forcing REF/ALT from genome"
  ./bin/plink2 --pfile "$STUDY" \
    --fa "$FA" --normalize --ref-from-fa force \
    --autosome --snps-only just-acgt --maf "$MAF" --geno "$GENO" \
    --set-all-var-ids 'chr@:#:$r:$a' --sort-vars \
    --make-pgen --threads "$THREADS" --out "$OUTDIR/study_norm"
  ./bin/plink2 --pfile "$OUTDIR/ref_pgen" \
    --fa "$FA" --normalize --ref-from-fa force \
    --autosome --snps-only just-acgt --maf "$MAF" --geno "$GENO" \
    --set-all-var-ids 'chr@:#:$r:$a' --sort-vars \
    --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_norm"
  STUDY="$OUTDIR/study_norm"; REF_PGEN="$OUTDIR/ref_norm"
else
  info "Standardizing IDs/QC (no FASTA provided)"
  ./bin/plink2 --pfile "$OUTDIR/ref_pgen" \
    --set-all-var-ids 'chr@:#:$r:$a' --sort-vars \
    --autosome --snps-only just-acgt --maf "$MAF" --geno "$GENO" \
    --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_named"
  ./bin/plink2 --pfile "$STUDY" \
    --set-all-var-ids 'chr@:#:$r:$a' --sort-vars \
    --autosome --snps-only just-acgt --maf "$MAF" --geno "$GENO" \
    --make-pgen --threads "$THREADS" --out "$OUTDIR/study_named"
  STUDY="$OUTDIR/study_named"; REF_PGEN="$OUTDIR/ref_named"
fi

# 2) Optional: drop ambiguous A/T & C/G
if [[ "$DROP_AMBIG" == "yes" ]]; then
  info "Dropping ambiguous palindromic SNPs"
  awk 'BEGIN{FS=OFS="\t"} $1!~/^#/{
    r=toupper($4); a=toupper($5);
    if(a ~ /,/) next;
    if(r!~/^[ACGT]$/ || a!~/^[ACGT]$/) next;
    if(!((r=="A"&&a=="T")||(r=="T"&&a=="A")||(r=="C"&&a=="G")||(r=="G"&&a=="C"))) print $3
  }' "$STUDY.pvar" > "$OUTDIR/keep_nonambig.ids"
  ./bin/plink2 --pfile "$STUDY"    --extract "$OUTDIR/keep_nonambig.ids" --make-pgen --threads "$THREADS" --out "$OUTDIR/study_noambig"
  ./bin/plink2 --pfile "$REF_PGEN" --extract "$OUTDIR/keep_nonambig.ids" --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_noambig"
  STUDY="$OUTDIR/study_noambig"; REF_PGEN="$OUTDIR/ref_noambig"
fi

# 3) Intersect IDs; align alleles; make shared BEDs
cut -f3 "$REF_PGEN.pvar" | sed '1d' | sort > "$OUTDIR/ref.ids"
./bin/plink2 --pfile "$STUDY" --extract "$OUTDIR/ref.ids" --sort-vars \
  --make-pgen --threads "$THREADS" --out "$OUTDIR/study_shared"
./bin/plink2 --pfile "$OUTDIR/study_shared" --make-bed --threads "$THREADS" --out "$OUTDIR/study_shared_b"
grep -v '^#' "$OUTDIR/study_shared_b.bim" | cut -f2 | sort > "$OUTDIR/shared.ids" || true
./bin/plink2 --pfile "$REF_PGEN" --extract "$OUTDIR/shared.ids" \
  --alt-allele force "$OUTDIR/study_shared_b.bim" 5 2 \
  --sort-vars --make-pgen --threads "$THREADS" --out "$OUTDIR/ref_shared"
./bin/plink2 --pfile "$OUTDIR/ref_shared" --make-bed --threads "$THREADS" --out "$OUTDIR/ref_shared_b"
NSHARED=$(wc -l < "$OUTDIR/shared.ids" || echo 0); info "Shared variants (pre-prune): $NSHARED"

# 3.5) Merge to get a "good by merge" SNP set
plink --bfile "$OUTDIR/ref_shared_b" --bmerge "$OUTDIR/study_shared_b" \
  --maf "$MAF" --geno "$GENO" --make-bed --out "$OUTDIR/merged_temp" --threads "$THREADS"
cut -f2 "$OUTDIR/merged_temp.bim" | sort > "$OUTDIR/merged_good.ids"

# 4) LD prune on REF only; apply same prune list to both
PRUNE_FLAGS=(--indep-pairwise "$PRUNE_WIN" "$PRUNE_STEP" "$PRUNE_R2")
[[ -n "$EXCL_RANGE" ]] && PRUNE_FLAGS+=(--exclude range "$EXCL_RANGE")

plink --bfile "$OUTDIR/ref_shared_b" "${PRUNE_FLAGS[@]}" \
  --extract "$OUTDIR/merged_good.ids" --threads "$THREADS" --out "$OUTDIR/$OUTPFX.prune"

plink --bfile "$OUTDIR/ref_shared_b"   --extract "$OUTDIR/$OUTPFX.prune.prune.in" --make-bed --threads "$THREADS" --out "$OUTDIR/ref_train"
plink --bfile "$OUTDIR/study_shared_b" --extract "$OUTDIR/$OUTPFX.prune.prune.in" --make-bed --threads "$THREADS" --out "$OUTDIR/study_project"

NPRUNE=$(wc -l < "$OUTDIR/$OUTPFX.prune.prune.in" || echo 0); info "Pruned markers (training set): $NPRUNE"

# 5) Train PCs on REF; project REF+STUDY on identical scale
./bin/plink2 --bfile "$OUTDIR/ref_train" --pca "$K" biallelic-var-wts \
  --out "$OUTDIR/$OUTPFX.ref.pcs" --threads "$THREADS"
./bin/plink2 --bfile "$OUTDIR/ref_train" --freq --out "$OUTDIR/$OUTPFX.ref"

PC_START=5; PC_END=$((K+4))  # eigenvec.var: CHROM ID MAJ NONMAJ PC1..PCK

./bin/plink2 --bfile "$OUTDIR/ref_train" \
  --read-freq "$OUTDIR/$OUTPFX.ref.afreq" \
  --score "$OUTDIR/$OUTPFX.ref.pcs.eigenvec.var" 2 4 header-read \
          variance-standardize list-variants \
  --score-col-nums ${PC_START}-${PC_END} \
  --out "$OUTDIR/$OUTPFX.ref.proj"

./bin/plink2 --bfile "$OUTDIR/study_project" \
  --read-freq "$OUTDIR/$OUTPFX.ref.afreq" \
  --score "$OUTDIR/$OUTPFX.ref.pcs.eigenvec.var" 2 4 header-read \
          variance-standardize list-variants \
  --score-col-nums ${PC_START}-${PC_END} \
  --out "$OUTDIR/$OUTPFX.study.proj"

# 6) Combine outputs for plotting (FID IID PC1..PCK on same scale)
awk -vK="$K" 'BEGIN{FS=OFS="\t"}
FNR==1{
  delete h
  for(i=1;i<=NF;i++) h[$i]=i
  if(!printed++){
    printf "#FID\tIID"; for(j=1;j<=K;j++) printf "\tPC%d", j; print ""
  }
  fid=(h["FID"]?h["FID"]:h["#FID"]); iid=h["IID"]
  pc=h["PC1_AVG"]?h["PC1_AVG"]:(h["PC1"]?h["PC1"]:0)
  if(!pc) for(i=1;i<=NF;i++) if($i~/^PC[0-9]+(_AVG)?$/){pc=i; break}
  next
}
{
  printf "%s\t%s", $(fid), $(iid)
  for(j=0;j<K;j++) printf "\t%s", $(pc+j)
  print ""
}' "$OUTDIR/$OUTPFX.ref.proj.sscore" "$OUTDIR/$OUTPFX.study.proj.sscore" \
> "$OUTDIR/$OUTPFX.combined.eigenvec"
cp "$OUTDIR/$OUTPFX.ref.pcs.eigenval" "$OUTDIR/$OUTPFX.combined.eigenval"

info "Done."
info "Outputs:"
info "  REF PCs:            $OUTDIR/$OUTPFX.ref.pcs.{eigenvec,eigenval,eigenvec.var}"
info "  REF allele freqs:   $OUTDIR/$OUTPFX.ref.afreq"
info "  REF projection:     $OUTDIR/$OUTPFX.ref.proj.sscore"
info "  STUDY projection:   $OUTDIR/$OUTPFX.study.proj.sscore"
info "  Combined PCs:       $OUTDIR/$OUTPFX.combined.{eigenvec,eigenval}"
