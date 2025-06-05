function  site_2_num = cd_site_inverse(matrix)  


site_a=matrix;
% 1 to site num
for i=1:size(site_a,2)
index=site_a(:,i)==1;
site_a(index,i)=i;
end

% matrix to num
site_2_num=[];
for i=1:size(site_a,1)
lin=sum(site_a(i,:));
site_2_num=[site_2_num;lin];
end

% delay site to num

k=size(site_a,2)+1;

index=site_2_num==0;
site_2_num(index)=k;

end