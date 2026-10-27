clc;
clear;
%Pridiction
xx=readmatrix('AlCoCrFeNi_FCC.xlsx','sheet','Sheet2','range','I2:V55211');

y=gh_ANN(xx);