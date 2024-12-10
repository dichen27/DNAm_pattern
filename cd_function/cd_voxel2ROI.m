function SST_264 =cd_voxel2ROI(SST_14,data_index,index_264)

data_index=[[1:length(data_index)]',data_index];
voxel_intersect=intersect(index_264(:,2),data_index(:,2));


SST_index=[];
for i=1:length(voxel_intersect)
    
    index=data_index(:,2)==voxel_intersect(i);
    
    lin=data_index(index,:);
    SST_index=[SST_index;lin];
end


mask_index=[];
for i=1:length(voxel_intersect)
    
index=index_264(:,2)==voxel_intersect(i);
lin=index_264(index,:);

  mask_index=[mask_index;lin];   
end



SST_14_264=[];
for i=1:size(SST_index,1)
    
lin=SST_14(:,SST_index(i,1));
SST_14_264=[SST_14_264,lin];
end


SST_264=[];
 for i=1:length(unique(index_264(:,1)))
     
 index=mask_index(:,1)==i; 
lin=nanmean(SST_14_264(:,index),2);

SST_264=[SST_264,lin];
 end
 
 
 
 
 
 
end










