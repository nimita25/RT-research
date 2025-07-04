ptid='9306087'; 
addpath('../pMBRT');
method=1;
gridDim = '111';
dr = 10;
ddr = 3;
px=2;nfrac=10;px0=px*nfrac;pvdr = 5;

% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
%load([folder ptid '/' ptid '.mat'],'cst','ct');
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);

x0=doseGrid.x;
y0=doseGrid.y;
z0=doseGrid.z;
cubeDim=doseGrid.dimensions;
x1=ct.x;
y=ct.y;
z=ct.z;
[x0,y0,z0] = ndgrid(x0,y0,z0);
[x1,y,z] = ndgrid(x1,y,z);

% 10,4
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\grid_113\res_9306087_2_1_2024-10-27-19-26.mat','d')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-29-08-03.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-11-04-00-52.mat','d')

% 10, 6
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-28-07-56.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-29-08-02.mat')

% 10, 3
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\grid_113\res_9306087_2_1_2024-11-01-08-12.mat','d')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-11-01-16-15.mat','d')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-11-03-21-31.mat','d')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-05-08-26.mat','d')


%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-13-02-04.mat','d','R')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-16-20-28.mat','d','R')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-16-21-02.mat','d','R')


%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-16-22-52.mat','d','R')

%load('output_9306087/res_9306087_8_4_2024-12-04-04-25.mat')
%load('output_9306087/res_9306087_8_2_2024-12-03-23-55.mat')
load('output_9306087/res_9306087_8_1_2024-12-03-21-18.mat')
%load('output_9306087/res_9306087_8_3_2024-12-04-00-58.mat')


% Load lattice info
% fname = strcat('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\LVertices\LVertices' ,gridDim, '_',num2str(dr),'_',num2str(ddr),'_9306087_',num2str(R(1)),num2str(R(2)),num2str(R(3)),'.mat');
%fname = strcat('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\LVertices111_',num2str(dr),'_',num2str(ddr),'_9306087_',num2str(R(1)),num2str(R(2)),num2str(R(3)),'.mat','Vvalley','Vpeak');
% load(fname);
load('9306087_lattice.mat');

d=d/(px*pvdr);
d=reshape(d,cubeDim);
d(d>1.1999)=1.1199;
ct.cubeHU{1}=interpn(x1,y,z,ct.cubeHU{1},x0,y0,z0);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;


tmp=cst([38],:);
tmp{1,4}{1}=[Vpeak;Vvalley];
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

for i=1:size(cst,1)
mask=zeros([nx ny nz]);
mask(cst{i,4}{1})=1;
mask=mask(idx,idy,idz);
cst{i,4}{1}=find(mask==1);
end
matRadGUI