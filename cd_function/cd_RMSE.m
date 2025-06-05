function rmse=cd_RMSE(y_real,y_prediect)
% Root Mean Squared Error，RMSE


%y_real = VDI_all(:,1);

%y_prediect = VDI_all(:,2);

% RMSE
rmse = sqrt(mean((y_real - y_prediect).^2)); 


end