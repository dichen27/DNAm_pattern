function [LD_cutoff,t_cutoff,voxel_index]=cd_LD_cutoff(R_2_moment,t_all,r_all,cutoff_d)

voxel_index_all=[1:length(t_all)]';


all=[R_2_moment,t_all,voxel_index_all];

d_all=cd_r2d(r_all);

index=abs(d_all)>cutoff_d;
all(index,:)=[];



LD_cutoff=all(:,1);
t_cutoff=all(:,2);
voxel_index=all(:,3);

Num_voxel=length(t_cutoff);





end