clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/ADNI/230211-cd_function/'));
working_dir='/Users/dichen/Documents/ADNI/';
cd(working_dir);

%% load ID_index
%ID_index=ADNIinformation;
%save('./ID_index.mat','ID_index');

load('./ID_index.mat')

ID_index_pan_str=cell2mat(table2array(ID_index(:,2)));
ID_index_ge_str=cell2mat(table2array(ID_index(:,3)));

%% logic ID_1919

ID_source=load('./ID_1919.mat');
ID_cell=ID_source.ID_1919;

ID_1919_logic=[];
for i=1:length(ID_cell)
i
ID_str=ID_cell{i};
pan=ID_str(1:12);
ge=ID_str(14:end);

index_matrix=[];
for j=1:size(ID_index_pan_str,1)
pan_lin=double(strcmp(pan, ID_index_pan_str(j,:)));
ge_lin=double(strcmp(ge, ID_index_ge_str(j,:)));
lin=[pan_lin,ge_lin];
index_matrix=[index_matrix;lin];
end
lin=[];
index=sum(index_matrix,2)==2;
if sum(index)==1
    lin=1;
else
    lin=0;
end


ID_1919_logic=[ID_1919_logic;[i,lin]];
end
save('./ID_1919_logic.mat','ID_1919_logic')


%% load ID_1919

ID_source=load('./ID_1905.mat');
ID_cell=ID_source.ID_1905;

ID_order_final=[];
for i=1:length(ID_cell)
i
ID_str=ID_cell{i};
pan=ID_str(1:12);
ge=ID_str(14:end);

index_matrix=[];
for j=1:size(ID_index_pan_str,1)
pan_lin=double(strcmp(pan, ID_index_pan_str(j,:)));
ge_lin=double(strcmp(ge, ID_index_ge_str(j,:)));
lin=[pan_lin,ge_lin];
index_matrix=[index_matrix;lin];
end

index=sum(index_matrix,2)==2;
ID_final=table2array(ID_index(index,1));

ID_order_final=[ID_order_final;ID_final];
end

%% load gender_self
%gender_self=ADNIcov;
%gender_self=table2array(gender_self)
%save('./gender_self.mat','gender_self');

load('./gender_self.mat')

gender_self_order=[];
for i=1:length(ID_order_final)

index=gender_self(:,1)==ID_order_final(i);
gender_lin=gender_self(index,2);

gender_self_order=[gender_self_order;[ID_order_final(i),gender_lin(1)]];
end

save('./gender_self_order.mat','gender_self_order')% 1:male 0:female
tabulate(gender_self_order(:,2))


%% load race_self
Race_self=xlsread('./Race_ADNI.xlsx');

Race_self_order=[];
for i=1:length(ID_order_final)

index=Race_self(:,1)==ID_order_final(i);
Race_lin=Race_self(index,3);

Race_self_order=[Race_self_order;[ID_order_final(i),Race_lin(1)]];
end

save('./Race_self_order.mat','Race_self_order')
tabulate(Race_self_order(:,end))

%% finished

close all; % close all figures
cd(working_dir);




%%











