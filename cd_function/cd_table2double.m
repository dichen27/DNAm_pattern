function ID_double= cd_table2double(ID_table)  

ID_array=table2array(ID_table);

ID_double=[];
for i=1:length(ID_array)
ID_double=[ID_double;str2double(cell2mat(ID_array(i)))];
end


end