%clc;clear;close all;
ptid = '9306087';
%folder=['D:\KUMC\codes 20240807\data\' ptid '\'];
%addpath(genpath('D:\KUMC\codes 20240807\matRad-master'));
%load([ptid '\' ptid '_c.mat'],'c');
%load([folder ptid '.mat'],'ct','cst');
folder = ['../' ptid '/'];
addpath(genpath('..\matRad-master'));
load([folder ptid '.mat'], 'ct', 'cst');
t_px = (1.2/1.1)*60; %total prescription dose
nfrac = 60; % number of fraction
px = t_px/nfrac; % prescription dose

%% Define target and OAR
ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy; Dmax<54Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy; Dmax<36Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy; Dmax<54Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy; Dmax<54Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc; Dmax<60Gy
ctv = cell(1,1);
ctv{1} = ctv1; %row indices of Dij corresponding to target/tumor
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);  %row indices of Dij corresponding to OAR
    n_oar(i) = numel(oar{i}); %number of voxels in each OAR
end
c = [ctv; {body}; oar;]; %cell containing row indices corresponding to OAR, tumor, body



x=ct.x;
y=ct.y;
z=ct.z;
[x,y,z] = ndgrid(x,y,z);


name1=['./output/BED-ADMM-2024-08-31-18-15.mat']; % ARC
name2=['./output/BED-ADMM-2024-08-31-18-50.mat']; % multi-IMPT
% name1=['./output/BED-ADMM-2024-08-31-18-50.mat']; % NC-BED
% name2=['./output/BED-ADMM-2024-08-31-19-23.mat']; % NC-Dose
% name3=['./output/BED-ADMM-2024-08-31-18-15.mat']; % Conv

% load(name3,'d');
% dconv = d;
% nn = 1:numel(dconv);
% for i = 1:numel(c)
%     if i ~= 2
%         nn = setdiff(nn,c{i});
%     end
% end


load(name2,'d');

% Get mean doses for further calculations
TD = zeros(size(d,1),1);
unique_fields = size(d,2);
for nf = 1:nfrac
    TD = TD+d(:,mod(nf-1,unique_fields)+1);
end
%TD = sum(d,2);
MD = TD/nfrac;
d = MD;
%d(nn) = dconv(nn);

d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,ct.cubeDim);
d(d>1.199)=1.199;
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([38],:);
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];%(!)

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

% for i=1:size(cst,1)
% mask=zeros([nx ny nz]);
% mask(cst{i,4}{1})=1;
% mask=mask(idx,idy,idz);
% cst{i,4}{1}=find(mask==1);
% end
matRadGUI
