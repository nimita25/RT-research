ptid='7119049'; 
%addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
gridDim = '113';
px = 2; % prescription dose
nfrac = 30; % number of fraction
mu_min = 5; % MMU threshold

% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
load([folder 'dij_' ptid '_doseGrid113.mat']);

%% Define target and OAR
ctv1=cst{9,4}{1};   % ptv60, 2*30
body=cst{1,4}{1};
N_oar=3; % (!)
oar=cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
oar(4)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};


ctv = cell(1,1);
ctv{1} = ctv1;
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);
    n_oar(i) = numel(oar{i});
end
c = [ctv; {body}; oar;];
%==============================================
N_c = numel(c);
n_c = zeros([N_c 1]);
for i = 1:N_c
    n_c(i) = numel(c{i});
end

for i=1:N_c
tmpCube=zeros(ct.cubeDim);
tmpCube(c{i}) = 1;
% interpolate cube
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
c{i} = VdoseGrid;
end


x0=doseGrid.x;
y0=doseGrid.y;
z0=doseGrid.z;
cubeDim=doseGrid.dimensions;
x=ct.x;
y=ct.y;
z=ct.z;
[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);


resolution=[1 1 3]; % mm
V = [cst{:,4}];
V = unique(vertcat(V{:}));
tmpCube    = zeros(ct.cubeDim);
tmpCube(V) = 1;
% interpolate cube
tmp_x=ct.resolution.x:resolution(1):ct.cubeDim(1)*ct.resolution.x;
tmp_y=(ct.resolution.y:resolution(2):ct.cubeDim(2)*ct.resolution.y)';
tmp_z=ct.resolution.z:resolution(3):ct.cubeDim(3)*ct.resolution.z;
DoseNumber = find(matRad_interp3((1:ct.cubeDim(1))*ct.resolution.x,(1:ct.cubeDim(2))*ct.resolution.y,(1:ct.cubeDim(3))*ct.resolution.z,...
    tmpCube,tmp_x,tmp_y,tmp_z,'nearest'));

% load('Results_7119049/res0602_7119049_3_shifts_1_RND_50_2025-06-26-13-30.mat','d')
load('Results_7119049/res0602_7119049_3_shifts_4_QC_50_2025-06-26-21-22.mat','d')
% load('Results_7119049/res0602_7119049_3_shifts_4_multiQC_50_2025-06-27-11-39.mat','d')

% d1 = zeros(prod(doseGrid.dimensions),1);
% d1(DoseNumber) = d(DoseNumber);
% d = d1;

d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,cubeDim(1), cubeDim(2), cubeDim(3));
d(d>1.199)=1.199;
% d=interpn(x0,y0,z0,d,x,y,z);
ct.cubeHU{1}=interpn(x,y,z,ct.cubeHU{1},x0,y0,z0);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([9],:);
tmp{1,4}{1}=c{1};
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];%(!)

nx=cubeDim(1);
ny=cubeDim(2);
nz=cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

ct.cubeDim=[numel(idx) numel(idy) numel(idz)];

ct.x=doseGrid.x(idx);
ct.y=doseGrid.y(idy);
ct.z=doseGrid.z(idz);
ct.cubeHU{1}=ct.cubeHU{1}(idx,idy,idz);
%resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,idz);

% for i=1:size(cst,1)
% mask=zeros([nx ny nz]);
% mask(cst{i,4}{1})=1;
% mask=mask(idx,idy,idz);
% cst{i,4}{1}=find(mask==1);
% end
matRadGUI