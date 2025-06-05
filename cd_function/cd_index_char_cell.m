function index=cd_index_char_cell(ID_ADHD,ID_select)
% size of bb:  94 * 94

index=[];
for i=1:length(ID_ADHD)
    
lin=strcmp(ID_ADHD{i},ID_select);

index=[index;lin];
end
index=logical(index);% double to logical
% sum(index)



end