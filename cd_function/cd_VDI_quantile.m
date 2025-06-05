function VDI_quantile=cd_VDI_quantile(VDI)


VDI_all=[[1:length(VDI)]',VDI];
standard_quantile=[1:length(VDI)]'/length(VDI);
VDI_sort=[sortrows(VDI_all,2),standard_quantile];% order,VDI, quantile


% re-sort to original
VDI_quantile_sort=sortrows(VDI_sort,1);
VDI_quantile=VDI_quantile_sort(:,end);

end