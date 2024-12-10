function data_zscore=cd_zscore(data_no_nan,std_want)

m=mean(data_no_nan);
s=std(data_no_nan);
xi=1/(std_want);


data_zscore=[];
for i=1:length(data_no_nan)
    
lin=(data_no_nan(i)-m)/(xi*s);

data_zscore=[data_zscore;lin];
end


end