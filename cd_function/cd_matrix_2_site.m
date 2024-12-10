function site=cd_matrix_2_site(site_matrix)


num=[1:size(site_matrix,2)];

site=[];
for i=1:size(site_matrix,1)
index=site_matrix(i,:)==1;

if sum(index)==0
lin=size(site_matrix,2)+1;
else
lin=num(index);
end

site=[site;lin];
end


end