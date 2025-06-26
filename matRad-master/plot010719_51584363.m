
clear

load 51584363.mat
x0=ct.x;
y0=ct.y;
z0=ct.z;

load 51584363_fine.mat
x=ct.x;
y=ct.y;
z=ct.z;

ctdim=ct.cubeDim;

[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);

load 51584363_dose_v2.mat dmean_photon dmean_proton dmean_hybrid
d=interpn(x0,y0,z0,dmean_proton/2,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([23;21;8;6;7;9;10],:);
cst=tmp;

ctv=cst{2,4}{1};
cst{3,4}{1}=setdiff(cst{3,4}{1},ctv);
cst{4,4}{1}=setdiff(cst{4,4}{1},ctv);
cst{5,4}{1}=setdiff(cst{5,4}{1},ctv);
cst{6,4}{1}=setdiff(cst{6,4}{1},ctv);
cst{7,4}{1}=setdiff(cst{7,4}{1},ctv);


% cst{1,5}.visibleColor=[1 0 0];
cst{1,5}.visibleColor=[1 0 1];
cst{2,5}.visibleColor=[1 1 1];
cst{3,5}.visibleColor=[0 0 0];
cst{4,5}.visibleColor=[0 0 1];
cst{5,5}.visibleColor=[0 0 1];
cst{6,5}.visibleColor=[0 1 0];
cst{7,5}.visibleColor=[0 1 0];
% cst{1,5}.visibleColor=[1 0 1];
% cst{2,5}.visibleColor=[1 1 1];
% cst{3,5}.visibleColor=[0 0 1];
% cst{4,5}.visibleColor=[0 0 1];
% cst{5,5}.visibleColor=[0 1 0];
% cst{6,5}.visibleColor=[0 1 0];

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

if 1 % axis
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

for i=1:size(cst,1)
mask=zeros([nx ny nz]);
mask(cst{i,4}{1})=1;
mask=mask(idx,idy,idz);
cst{i,4}{1}=find(mask==1);
end

% figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.2 1.6]*30);set(gca,'FontSize', 18);
% figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.2 1]*25);set(gca,'FontSize', 18);
% figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.2 1.8]*28);set(gca,'FontSize', 18);
figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.2 0.7]*35);set(gca,'FontSize', 18);
% return
matRadGUI
