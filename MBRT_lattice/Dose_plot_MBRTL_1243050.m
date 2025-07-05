ptid='1243050'; 
addpath('../pMBRT');
method=1;
gridDim = '111';
dr = 10;
ddr = 3;
px=1.8;nfrac=28;px0=px*nfrac;pvdr = 5;

% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
%load([folder ptid '/' ptid '.mat'],'cst','ct');
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);

x0=doseGrid.x;
y0=doseGrid.y;
z0=doseGrid.z;
%cubeDim=doseGrid.dimensions;
cubeDim=doseGrid.dimensions;
x1=ct.x;
y=ct.y;
z=ct.z;
[x0,y0,z0] = ndgrid(x0,y0,z0);
[x1,y,z] = ndgrid(x1,y,z);


%load('output_1243050/res_1243050_3_1_2025-07-04-15-48.mat')
%load('output_1243050/res_1243050_3_4_2025-07-04-17-01.mat')
%load('output_1243050/res_1243050_3_1_2025-07-05-08-08.mat')
%load('output_1243050/res_1243050_3_2_2025-07-05-08-21.mat')
%load('output_1243050/res_1243050_3_3_2025-07-05-08-34.mat')
load('output_1243050/res_1243050_3_4_2025-07-05-08-45.mat')

% Load lattice info
load('1243050_lattice.mat');
load('1243050_c.mat')

%d=d(1:prod(cubeDim));
d=d/(px*pvdr);
d=reshape(d,cubeDim);
d(d>1.1999)=1.1199;
ct.cubeHU{1}=interpn(x1,y,z,ct.cubeHU{1},x0,y0,z0);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;


N_oar = 4;
tmp = cell(N_oar+1,6);
count = 1;
for jj=[11 2 3 12 7]
for ii = 1:6
tmp{count,ii}=cst{jj,ii};
end
count=count+1;
end
tmp{1,4}{1} = [c{1};c{2}];
for ii = 2:N_oar+1
tmp{ii,4}{1}=[c{ii+2}];
end

cst=tmp;
cst{1,5}.visibleColor=[1 1 0];%(!)
cst{2,5}.visibleColor=[0 1 0];
cst{3,5}.visibleColor=[0.529 0.808 0.922];
cst{4,5}.visibleColor=[1 0.5 0];
cst{5,5}.visibleColor=[1 0 1];


% tmp=cst([11],:);
% tmp{1,4}{1}=[Vpeak;Vvalley];
% cst=tmp;
% cst{1,5}.visibleColor=[1 1 0];%(!)

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
resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,idz);

for i=1:size(cst,1)
mask=zeros([nx ny nz]);
mask(cst{i,4}{1})=1;
mask=mask(idx,idy,idz);
cst{i,4}{1}=find(mask==1);
end
matRadGUI