function data_dele=cd_dele_nan_row(data)

data(any(isnan(data), 2),:) = [];% dele any row NaN

data_dele=data;
end