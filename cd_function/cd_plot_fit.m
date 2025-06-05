
function cd_plot_fit(X_data,Y_data,cmd)
% g=gramm('x',X_data,'y',Y_data);
figure;
        plot(X_data,Y_data,'o','Color','k','MarkerSize',12);% oooo
 hold on;
        % Statistics (linear fit and plotting)
        b = [ones(length(X_data),1),X_data]\Y_data;
        
        
          x_fit=cmd(1,1:2);
%         x_fit = [min(X_data),max(X_data)];
%         x_fit=[0,24];
% %         plot(x_fit, x_fit * b(2) + b(1),'LineWidth',1.5);% line
          axis(cmd); % axis([xmin xmax ymin ymax]);
%         ylim([0,25]);% limit y only
end