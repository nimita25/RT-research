clear
ptid='7119049';
px=2;
folder=[ptid '\'];
addpath(genpath('C:\Users\nshinde\Desktop\matRad-master'));
load([ptid '_c.mat'],'c');
load([ptid '\dij_' ptid '_doseGrid113.mat']);
x0=doseGrid.x;
y0=doseGrid.y;
z0=doseGrid.z;
cubeDim=doseGrid.dimensions;

load([folder ptid '.mat'],'ct','cst');
x=ct.x;
y=ct.y;
z=ct.z;

[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);

name1=['res_' ptid '_3_3_3_ctc333']; % 3 mm
% name2=['res_' ptid '_3_3']; % convention
% name3=['res_' ptid '_3_4_2_4_0_0.1']; % uniform
name2=['res_' ptid '_3_17']; % convention
name3=['res_' ptid '_3_14_2_4_0_0.1_1']; % uniform
name4=['res_' ptid '_3_19']; % convention

load([folder name1 '.mat'],'d');
d=d/px;
d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
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

ct.x=x0(idx);
ct.y=y0(idy);
ct.z=z0(idz);
ct.cubeHU{1}=ct.cubeHU{1}(idx,idy,idz);
resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,idz);

for i=1:size(cst,1)
mask=zeros([nx ny nz]);
mask(cst{i,4}{1})=1;
mask=mask(idx,idy,idz);
cst{i,4}{1}=find(mask==1);
end

matRadGUI
