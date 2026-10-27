clc;
% Folder containing this script
projectDir = fileparts(mfilename('fullpath'));
% *************************************************************************
%  index of the atoms included in the database
% *************************************************************************
Li=3; Be=4; B=5; C=6; N=7; Na=11; Mg=12; Al=13; Si=14;
P=15; Ca=20; Sc=21; Ti=22; V=23; Cr=24; Mn=25; Fe=26;
Co=27; Ni=28; Cu=29; Zn=30; Ga=31; Ge=32; Sr=38; Y=39;
Zr=40; Nb=41; Mo=42; Tc=43; Ru=44; Rh=45; Pd=46; Ag=47;
Cd=48; In=49; Sn=50; Sb=51; La=57; Ce=58; Pr=59; Nd=60;
Pm=61; Sm=62; Eu=63; Gd=64; Tb=65; Dy=66; Ho=67; Er=68;
Tm=69; Yb=70; Lu=71; Hf=72; Ta=73; W=74; Re=75; Os=76;
Pt=78; Au=79; Pb=82; Bi=83;

data = readcell(fullfile(projectDir, 'AlCoCrFeNi_FCC.xlsx'),'sheet','sheet1','range','B2:D55211');
txt = [data(:,2), data(:,3)];
a_num = data(:,1);
num = cell2mat(a_num);
component=txt(:,1);% each componet in each alloy
fraction=txt(:,2);% each fraction in each alloy 
No=num;%number of elements in each alloy
Num=length(component);% number of types of alloys
Com=cell(Num,1);
features=zeros(Num,14);
vec = readmatrix(fullfile(projectDir, 'properties.xlsx'),'sheet','atomicsize','range','F2:F84');
vec_F = readmatrix(fullfile(projectDir, 'properties.xlsx'),'sheet','atomicsize','range','L2:L84');
for i=1:Num

% *************************************************************************
%  Calculate a series of alloy
% *************************************************************************

    Com(i) = regexp(component(i), '\s+', 'split');
    s=Com{i};
    Fra(i) = regexp(fraction(i), '\s+', 'split');
    s=Com{i};
    t=Fra{i};
    comp=zeros(No(i),1);
    ele_frac=zeros(No(i),1);
    for j=1:No(i)
        comp(j)=eval(cell2mat(s(j)));%transfer the data of components
        ele_frac(j)=eval(cell2mat(t(j)));%transfer the data of fractions
    end
    features(i,:)=calculatefeatures(comp, ele_frac);
    
end