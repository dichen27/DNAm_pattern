
clear;
clc;
close;

%addpath(genpath('/public/home/chendi/cd_function/'));
%working_dir='/public/home/chendi/Methylation_PRS/';
%cd(working_dir);

PRSice_dir='/home1/ISTBI_data/chendi/MDD_PRS/PRSice/';
cd(PRSice_dir)

%chmod -R 777 /home1/ISTBI_data/chendi/MDD_PRS/PRSice_v1.25/


%% MDD PRS

%code from Wang Xuefei
%code='R --file=PRSice_v1.25.R -q --args wd /home1/WangXF/mdd1 plink /home1/WangXF/PRSice_v1.25/plink_1.9_linux_160914 base /home1/WangXF/originGWAS/MDD2018_ex23andMe target /home1/WangXF/TARGET/imputed_plink_qc_merge_rs/imagen_merge_all_unique_filt_rs slower 0 sinc 0.005 supper 0.5 clump.snps T cleanup F no.regression T report.individual.scores T pheno.file /home1/WangXF/phenotext.txt'
    
wd_dir='/home1/ISTBI_data/chendi/MDD_PRS/result_MDD_PRS';% output
plink_dir='/home1/ISTBI_data/chendi/MDD_PRS/PRSice/plink_1.9_linux_160914';% plink
base_dir='/home1/ISTBI_data/chendi/MDD_PRS/MDD2018_ex23andMe';% GWAS result
target_dir='/home1/ISTBI_data/chendi/MDD_PRS/imputed_plink_qc_merge_rs/imagen_merge_all_unique_filt_rs';% imputed IMAGEN data
pheno_dir='/home1/ISTBI_data/chendi/MDD_PRS/PRSice/phenotype.txt';

cmd=['R --file=PRSice_v1.25.R -q --args wd ',wd_dir,' plink ',plink_dir,' base ',base_dir,' target ',target_dir,' slower 0 sinc 0.005 supper 0.5 clump.snps T cleanup F no.regression T report.individual.scores T pheno.file ',pheno_dir];


%cmd=['R --file=PRSice_v1.25.R -q --args wd ',wd_dir,' plink ',plink_dir,' base TOY_BASE_GWAS.assoc target TOY_TARGET_DATA slower 0 sinc 0.02 supper 0.3 covary F best.thresh.on.bar F report.individual.scores T']
unix(cmd);


%% finished




























