clc;
clear;
% Folder containing this script
projectDir = fileparts(mfilename('fullpath'));
%[num,txt] = xlsread('GH_FCCL12.xlsx','B2:D60000'); % read the componet and the fraction
data = readcell(fullfile(projectDir, 'AlCoCrFeNi_FCC.xlsx'),'sheet','sheet1','range','B2:D55211');
txt = [data(:,2), data(:,3)];
a_num = data(:,1);
num = cell2mat(a_num);
features = readmatrix('AlCoCrFeNi_FCC.xlsx','sheet','sheet1','range','I2:V55211');
component=txt(:,1);% each componet in each alloy
fraction=txt(:,2);% each fraction in each alloy 
No=num;%number of elements in each alloy
Num=length(component);% number of types of alloys
normalize = zeros(Num,14);
%按列normalization

for i=1:14
    FlattenedData = features(:,i)';
    MappedFlattened = mapminmax(FlattenedData,0,1);
    normalize(:,i) = MappedFlattened';
end