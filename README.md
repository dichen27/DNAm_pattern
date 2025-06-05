 dir of IMAGEN= ./Methylation_Running
 dir of PPMI=./PPMI_Running
 dir of ADNI=./ADNI_Running
 dir of Customized Functions=./cd_function

 The main code：./Methylation_Running
（1）-（14）：Please read the code in order
（Reply1-*）：Response after the first round of peer review
（Reply2-*）：Response after the second round of peer review


Result of 175 modules: ./Methylation_Running/Module175
Result of 18 Clusters: ./Methylation_Running/Cluster18

DNAm 18 Cluster pattern should be reordered as: 4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17

###########################################################################################
Catalogue of the IMAGEN 

(1)Methylation_regress.m
# The DNAm data were regressed out from the covariates for WVCNA analysis.

(2-1)Power_age14.R
# Annotation: Calculate the WGCNA parameter for DNAm at age 14, i.e., the power.

(2-2)Power_age14.R
# Calculate the WGCNA parameter for DNAm at age 19, i.e., the power.

(3)WVCNA_Consensus.R
# WVCNA and Consensus network analysis

(4)age19-age14.m
# Calculate the difference between DNAm data at age 19 and age 14, i.e., DNAm₁₉ − DNAm₁₄.

(5)Heatmap_and_Clustering.R
# Plot a heatmap and perform cluster-18 analysis.

(6-1)Functional_Enrichment_Mantel_test_Cluster18.R
# Functional enrichment analysis of the 18 clusters.

(6-2)Tissue_specific_Enrichment_Cluster18.R
# Tissue specific enrichment analysis of the 18 clusters.

(7)data_Cluster18.R
# Prepare the data for the 18 clusters for subsequent analyses.

(8)data_Module175.R
# Prepare the data for the 175 Modules for subsequent analyses.

(9)Freesurfer_data.m
# Load the MRI data

(10-pre)Prepare_CCA_data.m
# Prepare the MRI data for the CCA analyses

(10)Run_CCA_cluster10.R
# Run the CCA analyses of 10 brain-related clusters

(10-2)MRI_age19vsage14.m
# Run the MRI analyses: age-19 vs. age-14

(11)Behavior_DNAm_age19.m
# Behavior analyses at age 19

(12)Behavior_DNAm_age14.m
# Behavior analyses at age 14

(13)PRS_With_behavior.m
# PRS analyses

(14)Demographic Characteristics.m
# Calculation for the Demographic Characteristics

(Reply1-1-1)Clustering_direct_WGCNA.R
# Stability analysis: performing WGCNA directly without conducting Consensus analysis.

(Reply1-1-2)data_module190.m
(Reply1-1-3)Heatmap_module190.R
(Reply1-1-4)Heatmap_cluster18_450k_850k.R
# Stability analysis to reply the reviewer

(Reply1-2)Enrichment_Module_0.R
# Functional enrichment analysis of the Module_0.

(Reply1-3)Run_CCA_methy_MRI_Nobrain_cluster8.R
# Stability analysis:  the CCA analyses of 8 No-brain clusters

(Reply1-4)Module175-MDD.m
# Calculate the correlations between the 175 modules and depressive symptoms.

(Reply1-5)Z-test_consensusVSnoconsensus.R
# Z-test of the Consensus result and No-Consensus result.

(Reply1-6)Behavior_Controled_PRS_age14.m
# Behavioral analysis after controlling for PRS at age 14.

(Reply1-7)Behavior_Controled_PRS_age19.m
# Behavioral analysis after controlling for PRS at age 19.

(Reply1-8)Regress_cluster18_final_pattern.m
# The cluster-based data were regressed out from the covariates for pattern.

(Reply1-9)Heatmap_cluster18_age14_age19.R
# plot cluster-based Heatmap

(Reply1-10)FirstPC_module.R
# First PC of the module-based analysis.

(Reply1-11)Enrichment_stability.R
# Stability analysis of the Functional enrichment analysis.

(Reply1-12)Genomic_location_Cluster.R
# Plot genomic location of the 18 DNAm clusters.

###########################################################################################
Catalogue of the ADNI

(1)ADNI_clean_data.R
# Prepare the DNAm data of ADNI.

(2)ADNI_mapping_Mantel_test.R
# Plot 18 CLuster of ADNI.

(3)ADNI_grouping.m
# make the ADNI data into 3 subgroups

(4)ADNI_demographic.m
# Calculation for the Demographic Characteristics of ADNI.

###########################################################################################
Catalogue of the PPMI

(1)PPMI_mapping_cluster18_Mantel_test.R
# Plot 18 CLuster of PPMI.

(2)PPMI_ID_grouping.m
# make the PPMI data into 3 subgroups

(3)PPMI_Validation_depression.m
# Validate the association between DNAm and depression in the PPMI dataset.

(4)(4)PPMI_Characteric.m
# Calculation for the Demographic Characteristics of PPMI.

###########################################################################################
