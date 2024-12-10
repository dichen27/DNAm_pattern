function index=cd_index_char(ID_ADHD,ID_select)
% size of bb:  94 * 94

index=[];
for i=1:size(ID_ADHD,1)
    
lin=strcmp(ID_ADHD(i,:),ID_select);

index=[index;lin];
end
index=logical(index);% double to logical
% sum(index)



end