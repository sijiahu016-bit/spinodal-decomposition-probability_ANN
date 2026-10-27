function[features]=calculatefeatures(comp, ele_frac)
% *************************************************************************
%  Matlab Code to Compute the average radius and its deviation, the average
%  melt temperature and its deviation, the average mixing enthalpy and its
%  deviation, and the ideal entropy,the average electronal negativity and
%  its deviation, the average VEC and its deviation, the radio of
%  Scorr/Sid, the bulk modulus and its deviation.
% *************************************************************************
data1_ld = readmatrix('properties.xlsx','sheet','atomicsize','range','c2:H84');   %read the database
data2_ld = readmatrix('OQMD_Enthalpy.xlsx','sheet','sheet1','range','c3:CG85');  %read the mixing enthalpy in the database
Orthomatrix = readmatrix('Orthogonalization_matrix.xlsx','sheet','sheet1','range','a1:I9'); %read the orthomatrix
radius = data1_ld(:,1);  %the atomic size
tm = data1_ld(:,2);    %the melting point
elec_nega = data1_ld(:,3); %the Pauling electronic negativity
vec = data1_ld(:,4); %VEC
bulk = data1_ld(:,6);%bulk modulus
enthalpy = data2_ld;   %the mixing enthalpy
ntotal = length(comp); %the number of elements
ele_size=zeros(ntotal,1);
ele_temp=zeros(ntotal,1);
ele_elec_nega=zeros(ntotal,1);
VEC=zeros(ntotal,1);
for i=1:ntotal
    ele_size(i)=radius(comp(i));
    ele_temp(i)=tm(comp(i));
    ele_elec_nega(i)=elec_nega(comp(i));
    VEC(i)=vec(comp(i));
end  
% *************************************************************************
%   Calculate the fractions of each element 
% *************************************************************************
sum=0;
for i=1:ntotal
    sum=sum+ele_frac(i);
end
ele_frac=ele_frac/sum; 
% *************************************************************************
%   Calculate the mean redius and delta  
% *************************************************************************
r_mean=0;                              
for  i=1:ntotal
    r_mean=r_mean+ele_size(i)*ele_frac(i);      %the average atomic radius
end
delta=0;
for  i=1:ntotal
    delta=delta+ele_frac(i)*(1-ele_size(i)/r_mean)^2;
end
delta=sqrt(delta);  % the delta

% *************************************************************************
%   Calculate the average melting temperature and its deviation 
% *************************************************************************
TM=0;
for i=1:ntotal
    TM=TM+ele_frac(i)*ele_temp(i); % the average melting temperature
end
DTM=0;
for i=1:ntotal
    DTM=DTM+ele_frac(i)*(ele_temp(i)-TM)^2;
end
DTM=sqrt(DTM); % the deviation of melting temperature

% *************************************************************************
%   Calculate the average electronic negativity and its deviation 
% *************************************************************************
Mean_elecnega=0;
for i=1:ntotal
    Mean_elecnega=Mean_elecnega+ele_elec_nega(i)*ele_frac(i);
end
D_elecnega=0;
for i=1:ntotal
    D_elecnega=D_elecnega+ele_frac(i)*(ele_elec_nega(i)-Mean_elecnega)^2;
end
D_elecnega=sqrt(D_elecnega);

% *************************************************************************
%   Calculate the average electronic negativity and its deviation 
% *************************************************************************
MVEC=0;
for i = 1:ntotal
    MVEC=MVEC+ele_frac(i)*VEC(i);
end
D_VEC=0;
for i = 1:ntotal
    D_VEC=D_VEC+ele_frac(i)*(VEC(i)-MVEC)^2;
end
D_VEC=sqrt(D_VEC);

% *************************************************************************
%   Calculate the total mixing enthalpy 
% *************************************************************************
ME=0;
for i=1:ntotal-1
    for j=i+1:ntotal
        ME=ME+4*ele_frac(i)*ele_frac(j)*enthalpy(comp(i),comp(j));
    end
end   %the mixing enthalpy

% *************************************************************************
%   Calculate the deviation of mixing enthaply 
% *************************************************************************
DME=0;
for i=1:ntotal-1
    for j=i+1:ntotal
        DME=DME+ele_frac(i)*ele_frac(j)*(enthalpy(comp(i),comp(j))-ME)^2;
    end
end
DME=sqrt(DME); %the deviation of mixing enthaply

% *************************************************************************
%   Calculate the ideal mixing entropy
% *************************************************************************
Sid=0;
for i=1:ntotal
    Sid=Sid-ele_frac(i)*log(ele_frac(i));
end



% *************************************************************************
%   Calculate the Scorr/Sid
% *************************************************************************
%radio=ScorrSid(comp,ele_frac,TM,data1_ld,data2_ld);
 %atomic radius
 c=ele_frac;

    %atomic lattice constant of FCC structure
a=ele_size*2^1.5;

    %lattice constant
an=3.585;

    %coordination number of FCC
CN=12;

    %calculate w(i,j), A(i,j), the second term of equation 6
cjwA=0;
for  i=1:ntotal
    cjw(i)=0;
    cjA(i)=0;
    for j=1:ntotal
        x(i,j)=ele_size(i)/ele_size(j);
        w(i,j)=2*pi*(1-(ele_size(i)*(ele_size(i)+2*ele_size(j)))^0.5/(ele_size(i)+ele_size(j)));
        A(i,j)=2*pi*x(i,j)/(x(i,j)+1)^2/(x(i,j)*(x(i,j)+2))^0.5;
        cjw(i)=cjw(i)+c(j)*w(i,j);
        cjA(i)=cjA(i)+c(j)*A(i,j);
    end
    cjwA=cjwA+a(i)*c(i)*cjw(i)/cjA(i);
end  

aici=0;
temp1=0;
for i=1:ntotal
    aici=aici+a(i)*c(i);
    temp2=0;
    for j=1:ntotal
        temp2=temp2+c(j)*A(i,j);
    end
    temp1=temp1+a(i)*c(i)/temp2;
end
tem=(aici+cjwA-an)/temp1;%temp is 4*pi*n/N

    %calculate residual strain e(i)
for i=1:ntotal
    e(i)=(cjw(i)-tem)/cjA(i);
end

    %calculate ideal residual strain e(i)
n=0;
for i=1:ntotal
    for j=1:ntotal
        n=n+c(i)*c(j)*CN*(1-(x(i,j)*(x(i,j)+2))^0.5/(x(i,j)+1));
    end
end
n=n/2;
for i=1:ntotal
    ideal_e(i)=(cjw(i)-4*pi*n/CN)/cjA(i);
end

    %calculate root mean residual strain e(i)
root_e=0;
root_ie=0;
for i=1:ntotal
    root_e=root_e+c(i)*e(i)^2;
    root_ie=root_ie+c(i)*ideal_e(i)^2;
end
root_e=root_e^0.5;
root_ie=root_ie^0.5;

    %Lattice constant
a_d=0;
r_ave=0;
for i=1:ntotal
    r_ave=r_ave+c(i)*ele_size(i)*(1+ideal_e(i));
    a_d=a_d+c(i)*ele_size(i)*2^1.5*(1+ideal_e(i));
end

    %volume per atom
V=4/3*pi*(r_ave)^3*1e-30;

    %average bulk modulus B
B_ave=0;
B=zeros(ntotal,1);
for i=1:ntotal
    B(i)=data1_ld(comp(i),6);
end
for i=1:ntotal
    
    B_ave=B_ave+c(i)*B(i);
end
D_Bulk=0;
for i = 1:ntotal
    D_Bulk=D_Bulk+c(i)*(B(i)-B_ave)^2;
end
D_Bulk=sqrt(D_Bulk);
B_ave=B_ave*1e9;



    %Boltzmann constant
k=1.380645*1e-23;

    %elastic energy fluctuation delta elastic energy/kT
x_e=3*2^0.5*root_ie*(B_ave*V/(k*TM))^0.5;




    %mixing enthalpy
H_ave=0;
t=0;
for i=1:ntotal
    for j=1:ntotal
        if j~=i
            H_ave=H_ave+enthalpy(comp(i),comp(j));
        end
    end
end
H_ave=H_ave/(ntotal^2-ntotal);

D_H=0;
for i=1:ntotal
    for j=1:ntotal  
       if j~=i
            D_H=D_H+(data2_ld(comp(i),comp(j))-H_ave)^2; 
           end
    end
end
D_H=(D_H/(ntotal^2-ntotal))^0.5;

temp=0;
for i=1:ntotal
    for j=1:ntotal  
       if j~=i
            temp=temp+c(i)*c(j)*(data2_ld(comp(i),comp(j))-ME)^2; 
           end
    end
end

    %avogadro's constant
NA=6.0221*1e23;

    %chemical bond fluctuation 
x_cb=2*(temp^0.5*1e3/NA/k/TM)^0.5;


    %ideal mixing entropy
S_id=0;
for i=1:ntotal
    S_id=S_id+c(i)*log(c(i));
end
S_id=-S_id;
Scorr=S_id+SE(x_e+x_cb);
ratio=Scorr/S_id;

% *************************************************************************
%   Calculate the eigenvalues of hessian matrix
% *************************************************************************
% gas constant R J/mol/K
R=8.3144598;

T1=TM*0.001; % Divide by 1000 to complete the unit conversion from J to KJ
%calculate Hessian matrix

Hessian=zeros(ntotal-1);
Hmix=zeros(ntotal);
for i=1:ntotal-1
    Hmix(i,ntotal)=enthalpy(comp(i),comp(ntotal));
    Hessian(i,i)=(c(i)+c(ntotal))*R*T1/(c(i)*c(ntotal))-8*Hmix(i,ntotal);
    for j=i+1:ntotal-1
        Hmix(i,j)=enthalpy(comp(i),comp(j));
        Hmix(j,ntotal)=enthalpy(comp(j),comp(ntotal));
        if Hmix(i,j)==0
           Hmix(i,j)=0;
        end
        Hmix(j,i)=Hmix(i,j);
        Hmix(ntotal,i)=Hmix(i,ntotal);
        Hmix(ntotal,j)=Hmix(j,ntotal);
        Hessian(i,j)=4*(Hmix(i,j)-Hmix(i,ntotal)-Hmix(j,ntotal))+R*T1/c(ntotal);
        Hessian(j,i)=Hessian(i,j);
    end
end
%Matrix arithmetic
Ortho_Hessian=Orthomatrix(1:(ntotal-1),1:(ntotal-1))'*Hessian*Orthomatrix(1:(ntotal-1),1:(ntotal-1));
d=det(Ortho_Hessian);
[V,D]=eig(Ortho_Hessian);
B=Orthomatrix(1:(ntotal-1),1:(ntotal-1))*V;
A=diag(D);
Eig=min(A);

   

% *************************************************************************
%   Print all the features for machine learning
% *************************************************************************

features=zeros(1,14);
features(1)=r_mean; 
features(2)=delta;
features(3)=TM;
features(4)=DTM;
features(5)=ME;
features(6)=DME;
features(7)=Sid;
features(8)=Mean_elecnega;
features(9)=D_elecnega;
features(10)=MVEC;
features(11)=D_VEC;
features(12)=B_ave;
features(13)=D_Bulk;
features(14)=Eig;
