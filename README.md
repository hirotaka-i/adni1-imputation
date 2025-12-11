# ADNI1 genetic QC and population splitting

Install plink2 to ./bin folder `
```
mkdir -p ./bin
wget -O ./bin/plink2.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip && unzip -o ./bin/plink2.zip -d ./bin && chmod +x ./bin/plink2
```
Download chain file for liftOver
```
wget -O data/hg18ToHg38.over.chain.gz http://hgdownload.soe.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
````
Install modules to run the following scripts on cluster
```
module load ucsc #for liftOver
module load plink/1.9
```
## Required references
* hg38.fa.gz (fasta file) from [UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
* AJ reference plink binary from [GSE23636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23636)
* 1000 genome phase 3 plink binary + population labels filtered to biallelic snps on autosomes with a MAF > 0.01, geno > 0.95 and hwe > 1e-6. Also, pallindromes and long LD regions were excluded. ([Processing steps](https://github.com/hirotaka-i/1kg_ref/blob/main/main.ipynb))

## EUR
Start with removing the noisy variants from the original ADNI1 hg18 data. 
Then, align to hg38, standardize the variants IDs to chr:pos:ref:alt
```
awk '$5 == "0" || $6 == "0" || $5 == $6 {print $2}' \
  ../resources/ADNI/ADNI1_hg19/ADNI_cluster_01_forward_757LONI.bim \
  > temp/bad_snps.txt

./bin/plink2 \
  --bfile ../resources/ADNI/ADNI1_hg19/ADNI_cluster_01_forward_757LONI \
  --exclude temp/bad_snps.txt \
  --make-bed \
  --out temp/ADNI_start

bash code/01_prepare_variants.sh \
        --bfile temp/ADNI_start \
        --fa ../resources/liftover_ref/hg38.fa.gz \
        --chain data/hg18ToHg38.over.chain.gz \
        --outdir temp/prep_vars_out \
        --threads 2
```
sample QCs (call rate, sex, heterozygosity, relatives)
```
bash code/02_qc_samples.sh \
        --pfile temp/prep_vars_out/hg38_prepped \
        --outdir temp/qc_samples_out \
        --threads 2
```
Population splitting (Map study samples to 1kg reference panel PC projection)
```
bash code/03_pca_with_ref_score.sh \
        --study-pfile temp/qc_samples_out/final_keep \
        --ref-bfile ../resources/1kg_p3/all_hg38_filtered_chrpos \
        --fa ../resources/liftover_ref/hg38.fa.gz \
        --drop-ambig yes \
        --outdir temp/merge_ref_proj_out \
        --threads 2
python code/03b_plot_pop_and_split.py \
        --pc-prefix temp/merge_ref_proj_out/study_vs_ref.combined \
        --ref-label ../resources/1kg_p3/all_hg38_filtered_chrpos_pop.txt \
        --ref-label-col Population \
        --split-method mahalanobis \
        --out-prefix temp/pop_split_out/with_1kg_mah
```
```
InfPop
REF      2573
EUR       621
AAC        30
AMR        21
OTHER       9
SAS         1
```

Output figures are in `temp/pop_split_out`
![Whole population split](temp/pop_split_out/with_1kg_mah_Whole_group_PC1xPC2.png)
![EUR population split](temp/pop_split_out/with_1kg_mah_EUR_PC1xPC2.png)

# Separate EUR and prep for imputation
```
bash code/04_qc_split.sh \
  --pfile temp/prep_vars_out/hg38_prepped \
  --keep-samples temp/pop_split_out/with_1kg_mah_EUR.list \
  --fa ../resources/liftover_ref/hg38.fa.gz \
  --snps-only yes \
  --vcf-out yes \
  --outdir temp/qc_EUR_splitted \
  --geno-thres 0.1 \
  --maf-thres 0.005 \
  --vcf-out yes \
  --threads 2
```

Archive the temp/qc_EUR_splitted/vcf/chr*.vcf.gz output files for download
```
cd temp/qc_EUR_splitted/vcf && tar -czvf ../preimpute_vcf.tar.gz chr*.vcf.gz 
```

--> Michigan Imputation server upload with reference panel: 1000 Genomes Phase 3 v5 (GRCh38/hg38) x30