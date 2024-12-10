function [phenotype_no_nan,num_nan]=cd_nan(phenotype)

num_nan=sum(isnan(phenotype));
phenotype(isnan(phenotype))=nanmean(phenotype);


phenotype_no_nan=phenotype;



end