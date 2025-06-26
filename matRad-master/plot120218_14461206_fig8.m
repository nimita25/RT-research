
clear

ptid='14461206';
folder=[pwd '\' ptid '\'];

load([folder ptid '.mat'],'ct');
x0=ct.x;
y0=ct.y;
z0=ct.z;

load([folder ptid '_fine.mat'],'ct','cst');
x=ct.x;
y=ct.y;
z=ct.z;

ctdim=ct.cubeDim;

[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);

if 0
load([folder 'Gmin5_mu0001_m32_n200_wmax100_d.mat'],'d');
tmp=cst([14;15;17],:);
cst=tmp;
else
%     load([folder 'robust_Gmin5_mu0001_n200_v2_d.mat'],'d');
    load([folder 'robust_Gmin5_mu0001_m20_n200_v2_d.mat'],'d');
%     load([folder 'robust_Gmin5_mu0001_m33_n200_wmax100_d.mat'],'d');
d=d{1};
tmp=cst([9;15;17],:);
cst=tmp;
end

d=interpn(x0,y0,z0,d,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d/2;


cst{2,4}{1}=setdiff(cst{2,4}{1},cst{1,4}{1});
cst{3,4}{1}=setdiff(cst{3,4}{1},cst{1,4}{1});

cst{1,5}.visibleColor=[1 1 0];
cst{2,5}.visibleColor=[0 0 1];
cst{3,5}.visibleColor=[0 1 0];
% cst{4,5}.visibleColor=[0 1 0];
% cst{5,5}.visibleColor=[0 1 0];

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

if 0 % axis
idx=65:405;
idy=40:380;
idz=1:nz;
end

if 0 % sag
idx=100:340;
idy=1:ny;
idz=99:339;
end

ct.cubeDim=[numel(idx) numel(idy) numel(idz)];

ct.x=ct.x(idx);
ct.y=ct.y(idy);
ct.z=ct.z(idz);
ct.cubeHU{1}=ct.cubeHU{1}(idx,idy,idz);
resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,idz);

% figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.0 1.2]*60);set(gca,'FontSize', 18);
matRadGUI
