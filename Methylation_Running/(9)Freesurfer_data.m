clear;
clc;
close;

working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);

%% load data
DataDir='./Freesurfer_IMAGEN_data/BL/aparc_GrayVol.txt';
Data_all=importdata(DataDir);
GrayVol_14=Data_all.data;


DataDir='./Freesurfer_IMAGEN_data/FU2/aparc_GrayVol.txt';
Data_all=importdata(DataDir);
GrayVol_19=Data_all.data;


DataDir='./Freesurfer_IMAGEN_data/BL/aparc_SurfArea.txt';
Data_all=importdata(DataDir);
SurfArea_14=Data_all.data;


DataDir='./Freesurfer_IMAGEN_data/FU2/aparc_SurfArea.txt';
Data_all=importdata(DataDir);
SurfArea_19=Data_all.data;


DataDir='./Freesurfer_IMAGEN_data/BL/aparc_ThickAvg.txt';
Data_all=importdata(DataDir);
ThickAvg_14=Data_all.data;


DataDir='./Freesurfer_IMAGEN_data/FU2/aparc_ThickAvg.txt';
Data_all=importdata(DataDir);
ThickAvg_19=Data_all.data;


DataDir='./Freesurfer_IMAGEN_data/BL/aseg_Volume.txt';
Data_all=importdata(DataDir);
SubcortexVol_14=Data_all.data;
ICV_14=[SubcortexVol_14(:,1),SubcortexVol_14(:,end)];

DataDir='./Freesurfer_IMAGEN_data/FU2/aseg_Volume.txt';
Data_all=importdata(DataDir);
SubcortexVol_19=Data_all.data;
ICV_19=[SubcortexVol_19(:,1),SubcortexVol_19(:,end)];

%% save
save([working_dir,'(9)Freesurfer_data.mat'],'GrayVol_14','GrayVol_19','SurfArea_14','SurfArea_19','ThickAvg_14','ThickAvg_19','SubcortexVol_14','SubcortexVol_19','ICV_14','ICV_19');



%% finished
cd(working_dir);






























