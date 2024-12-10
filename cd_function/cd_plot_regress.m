

function [x_regressed,y_regressed]=cd_plot_regress(X_data,Y_data,cov_plot)

cov_plot(:,find(sum(cov_plot)==0))=[];% dele the all 0 in cov

% X_data=data_asymmetry_score_l(:,1)
% Y_data=data_asymmetry_score_l(:,2)
% cov_plot=data_asymmetry_score_l(:,3:size(data_asymmetry_score_l,2))


cov_plot_one=[ones(size(cov_plot,1),1),cov_plot];


[b,bint,r,rint,stats] = regress(X_data,cov_plot_one);
x_regressed=X_data-cov_plot_one*b;%r



[b,bint,r,rint,stats] = regress(Y_data,cov_plot_one);
y_regressed=Y_data-cov_plot_one*b;%r


cd_plot(x_regressed,y_regressed);


% [r_plot,p_plot]=corr(op_x,op_y)
end