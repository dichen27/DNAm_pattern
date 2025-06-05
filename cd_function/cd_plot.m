
function cd_plot(X_data,Y_data)


[r,p]=corr(X_data,Y_data)
% g=gramm('x',X_data,'y',Y_data);
figure;
        plot(X_data,Y_data,'.','Color','k');% oooo
 hold on;
        % Statistics (linear fit and plotting)
        b = [ones(length(X_data),1),X_data]\Y_data;
        
   
        x_fit = [min(X_data),max(X_data)];
%         x_fit=[0,24];
        plot(x_fit, x_fit * b(2) + b(1),'LineWidth',1.5);% line
%         set(handles,'xtick',0:1:25) % handles����ָ������������ľ��
%         set(handles,'ytick',0:1:25) % handles����ָ������������ľ��
end