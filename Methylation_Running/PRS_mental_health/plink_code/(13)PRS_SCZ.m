
clear;
clc;
close;


%% MDD PRS
% fixed
PRSice_dir='/home1/ISTBI_data/chendi/MDD_PRS/PRSice/';
cd(PRSice_dir)
plink_dir='/home1/ISTBI_data/chendi/MDD_PRS/PRSice/plink_1.9_linux_160914';% plink
target_dir='/home1/ISTBI_data/chendi/MDD_PRS/imputed_plink_qc_merge_rs/imagen_merge_all_unique_filt_rs';% imputed IMAGEN data
pheno_dir='/home1/ISTBI_data/chendi/MDD_PRS/PRSice/phenotype.txt'; %fake phenotype

% Variable
wd_dir='/home1/ISTBI_data/chendi/SCZ_PRS/result';% output
base_dir='/home1/ISTBI_data/chendi/SCZ_PRS/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv';% GWAS result


cmd=['R --file=PRSice_v1.25.R -q --args wd ',wd_dir,' plink ',plink_dir,' base ',base_dir,' target ',target_dir,' slower 0 sinc 0.005 supper 0.5 clump.snps T cleanup F no.regression T report.individual.scores T pheno.file ',pheno_dir];
unix(cmd);


%% finished




























