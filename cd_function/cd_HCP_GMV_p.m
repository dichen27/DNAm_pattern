function p=cd_HCP_GMV_p(H2_real)

simulation_final_dir='/share/inspurStorage/home1/chendi/Documents/LD_regression/HCP_T1/';
load([simulation_final_dir,'(4)HCP_simulation_noise.mat']);



 H2=output(:,2);
 H2_sort=sortrows(H2,-1);
 
 p=[];
 for i=1:length(H2_real)
 p_lin=(length(find(H2_sort>=H2_real(i))))/length(H2);
 
p=[p; p_lin];
 end
 
 
 
end

