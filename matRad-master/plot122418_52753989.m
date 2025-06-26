
clear

ptid='52753989';
folder=[pwd '\' ptid '\'];

load([folder ptid '.mat'],'ct');
x0=ct.x;
y0=ct.y;
z0=ct.z;
ctdim0=ct.cubeDim;

load([folder ptid '_fine.mat'],'ct','cst');
x=ct.x;
y=ct.y;
z=ct.z;

ctdim=ct.cubeDim;

[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);

if 1
load([folder ptid '_res_hybridpb_0_d.mat'],'d1','d2');
% load([folder ptid '_res_hybridpb_10_d.mat'],'d1','d2');
% load([folder ptid '_res_hybridpb_30_d.mat'],'d1','d2');
% tmp=cst([7;9;11;17],:);
% cst=tmp;
cst=cst([22;28;27;24],:);
d=zeros(ctdim0);d=d(:);
for i=1:9
d=d+d1{i}+d2;
% d=d+d1{i};
% d=d+d2;
end
d=d/9;
d=reshape(d,ctdim0);
else
end
px=45;



ctv=cst{1,4}{1};
cst{2,4}{1}=setdiff(cst{2,4}{1},ctv);
cst{3,4}{1}=setdiff(cst{3,4}{1},ctv);
cst{4,4}{1}=setdiff(cst{4,4}{1},ctv);

d=interpn(x0,y0,z0,d,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d/px;


cst{1,5}.visibleColor=[1 1 0];
cst{2,5}.visibleColor=[0 0 0];
cst{3,5}.visibleColor=[0 0 1];
cst{4,5}.visibleColor=[0 1 0];

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

if 0 % axis
idx=60:300;
idy=70:250;
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

% figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.0 1.2]*70);set(gca,'FontSize', 18);
matRadGUI
