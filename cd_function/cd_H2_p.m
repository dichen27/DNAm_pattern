function p=cd_H2_p(H2_real,k)

simulation_final_dir='/home1/chendi/Documents/LD_regression/Simulation_final/';
load([simulation_final_dir,'(5-0)H2_noise_Sum.mat']);


 H2=H2_all(:,k);
 H2_sort=sortrows(H2,-1);
  p=(length(find(H2_sort>=H2_real)))/length(H2);






end

