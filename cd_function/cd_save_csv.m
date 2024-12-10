function cd_save_csv(save_dir,data)

% save_dir='/home1/chendi/Documents/CWAS_ABIDE/plot/FC_regressed_all.csv';

mydata=cd_num2table(data);
writetable(mydata,save_dir)
end