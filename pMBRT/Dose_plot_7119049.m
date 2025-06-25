ptid='7119049'; 
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=1;
gridDim = '113';
px=2;%nfrac=10;px0=px*nfrac;pvdr = 5;

% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde.KUMC\Desktop\pMBRT\';
load([ptid '_c.mat'],'c');
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);

x0=doseGrid.x;
y0=doseGrid.y;
z0=doseGrid.z;
cubeDim=doseGrid.dimensions;
x=ct.x;
y=ct.y;
z=ct.z;
[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);

%load('C:\Users\nshinde\Desktop\pMBRT\7119049\res_7119049_3_3_2_ctc357.mat','d')
load('C:\Users\nshinde.KUMC\Desktop\pMBRT\7119049\res_7119049_3_3_2_ctc355.mat','d')

d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,cubeDim);
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
resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,idz);

% for i=1:size(cst,1)
% mask=zeros([nx ny nz]);
% mask(cst{i,4}{1})=1;
% mask=mask(idx,idy,idz);
% cst{i,4}{1}=find(mask==1);
% end
matRadGUI