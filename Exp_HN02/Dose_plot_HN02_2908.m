%clc;clear;close all;
ptid = 'HN02';
%folder=['D:\KUMC\codes 20240807\data\' ptid '\'];
%addpath(genpath('D:\KUMC\codes 20240807\matRad-master'));
%load([ptid '\' ptid '_c.mat'],'c');
%load([folder ptid '.mat'],'ct','cst');
folder = ['../' ptid '/'];
addpath(genpath('C:\Users\nshinde\Desktop\matRad-master'));
load([folder ptid '.mat'], 'ct', 'cst');
t_px = 70; %total prescription dose
nfrac = 35; % number of fraction
px = t_px/nfrac; % prescription dose

%% Define target and OAR
ctv1 = cst{15,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{9,4}{1};
N_oar = 4;
oar = cell(N_oar,1);
oar(1) = cst{2,4}(1); % R Parotid V50<30Gy
oar(2) = cst{11,4}(1); % OralCavity Dmean<40Gy
oar(3) = cst{17,4}(1); % Oropharynx Dmax<20Gy
oar(4) = cst{16,4}(1); % Larynx Dmax<20Gy
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


name1=['./output/BED-ADMM-2024-10-17-11-00.mat']; % ARC
name2=['./output/BED-ADMM-2024-10-17-10-47.mat']; % multi-IMPT
% name1=['./output/BED-ADMM-2024-09-05-14-20.mat']; % NC-BED
% name2=['./output/BED-ADMM-2024-09-05-14-27.mat']; % NC-Dose
% name3=['./output/BED-ADMM-2024-09-05-14-10.mat']; % Conv

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

tmp=cst([15],:);
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];%(!)
% Create structure as needed by matRadGUI
tmp1.type = 'square deviation';
tmp1.penalty = 800;
tmp1.dose = 30;
tmp1.EUD = NaN;
tmp1.volume = NaN;
tmp1.robustness = 'none';
cst{1,6} = tmp1;

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

for i=1:size(cst,1)
mask=zeros([nx ny nz]);
mask(cst{i,4}{1})=1;
mask=mask(idx,idy,idz);
cst{i,4}{1}=find(mask==1);
end
matRadGUI
