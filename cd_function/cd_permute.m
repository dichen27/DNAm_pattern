function [phenotype_permute,permute_index]=cd_permute(phenotype)
% size of bb:  94 * 94


data_randi=[phenotype,randperm(length(phenotype))'];
permute_index=data_randi(:,2);
data_sort=sortrows(data_randi,size(data_randi,2));
phenotype_permute=data_sort(:,1);

end