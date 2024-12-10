function [mse]=cd_mse(data)

sample_size=length(data);
std_data=std(data);
mse=std_data/sqrt(sample_size);

end