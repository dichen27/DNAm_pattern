clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/LD_regression_Reply1/230211-cd_function/'));
working_dir='/Users/dichen/Documents/PPMI/preprocess/';
cd(working_dir);

%% load ID_index

%ID_index=PPMIMethn524forLONI030718;
%save('./ID_index.mat','ID_index');

load('./ID_index.mat')

ID_index_pan_str=cell2mat(table2array(ID_index(:,2)));
ID_index_ge_str=cell2mat(table2array(ID_index(:,3)));



%% load ID_524

ID_source=load('./ID_524.mat');
ID_cell=ID_source.ID_524;

ID_order_final=[];
for i=1:524

ID_str=ID_cell{i};
pan=ID_str(1:12);
ge=ID_str(14:end);

index_matrix=[];
for j=1:524
pan_lin=double(strcmp(pan, ID_index_pan_str(j,:)));
ge_lin=double(strcmp(ge, ID_index_ge_str(j,:)));
lin=[pan_lin,ge_lin];
index_matrix=[index_matrix;lin];
end

index=sum(index_matrix,2)==2;

ID_final=cd_table2double(ID_index(index,1));

ID_order_final=[ID_order_final;ID_final];
end

%% load gender_self

gender_self=xlsread('./Gender.xlsx');

gender_self_order=[];
for i=1:524

index=gender_self(:,1)==ID_order_final(i);
gender_lin=gender_self(index,2)+1;

gender_self_order=[gender_self_order;[ID_order_final(i),gender_lin]];
end

save('./gender_self_order.mat','gender_self_order')

%% load race_self

Race_self=xlsread('./Race_PPMI.xlsx');

Race_self_order=[];
for i=1:524

index=Race_self(:,1)==ID_order_final(i);
Race_lin=Race_self(index,2);

Race_self_order=[Race_self_order;[ID_order_final(i),Race_lin]];
end

save('./Race_self_order.mat','Race_self_order')

%% finished

close all; % close all figures
cd(working_dir);




%%











