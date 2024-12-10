function cd_save_R(save_dir,data)


% save_dir='/home1/chendi/Documents/CWAS_ABIDE/plot/FC_regressed_all.rds';

test = fopen(save_dir,'wb');
fwrite(test,data,'double');
fclose(test);

end