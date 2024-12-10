function cd_plot_all(X_data,Y_data,cmd)

X_Y_num=[X_data,Y_data,ones(length(X_data),1)];

%% calculate the repeat num

for i=1:length(X_data)
for j=i+1:length(X_data)
      if  X_Y_num(j,1)==X_Y_num(i,1) &  X_Y_num(j,2)==X_Y_num(i,2)
      X_Y_num(i,3)=X_Y_num(i,3)+1;
      end
end
end

% if 14+(max(X_Y_num(:,3))-1)*3>30  %if MarkerSize too large: error!!!
%     error('cd27 oooo is out !!!')
% end

%%
% g=gramm('x',X_data,'y',Y_data);% main function
figure;

X_Y_num_reduce=X_Y_num;
for i=1:max(X_Y_num(:,3))
    
% main plot
 plot(X_Y_num_reduce(:,1),X_Y_num_reduce(:,2),'o','Color','k','MarkerSize',14+3*(i-1));% oooo
 
% reduce
 for j=1:size(X_Y_num_reduce,1)
    X_Y_num_reduce(j,3)=X_Y_num_reduce(j,3)-1; 
 end
     index=X_Y_num_reduce(:,3)==0;
      X_Y_num_reduce(index,:)=[];     
      
hold on;
end
   
       
b = [ones(length(X_data),1),X_data]\Y_data; % Statistics (linear fit and plotting)
           
% x_fit = [min(X_data),max(X_data)];
% x_fit=[0,25];
x_fit=cmd(1,1:2);
plot(x_fit, x_fit * b(2) + b(1),'LineWidth',1.5);% line
   axis(cmd); % axis([xmin xmax ymin ymax]);
%        ylim([0,25]);% limit y only
end