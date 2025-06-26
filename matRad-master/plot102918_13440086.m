
clear

load 13440086.mat
x0=ct.x;
y0=ct.y;
z0=ct.z;

ptv=cst{12,4}{1};
load 13440086_dose.mat d_robust_hybrid d_robust_hybrid_proton d_robust_hybrid_photon
md=mean(d_robust_hybrid{1}(ptv));
md1=mean(d_robust_hybrid_proton{1}(ptv));
md2=mean(d_robust_hybrid_photon{1}(ptv));

r1=30*md1/(md1+md2)
r2=30*md2/(md1+md2)

load 13440086_fine.mat
x=ct.x;
y=ct.y;
z=ct.z;

ctdim=ct.cubeDim;

[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);

load 13440086_dose.mat d_robust_hybrid d_robust_hybrid_proton d_robust_hybrid_photon
% md=mean(d_robust_hybrid{1}(ptv));
% md1=mean(d_robust_hybrid_proton{1}(ptv));
% md2=mean(d_robust_hybrid_photon{1}(ptv));

d=interpn(x0,y0,z0,d_robust_hybrid_photon{6},x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([12;7;6;2],:);
cst=tmp;

% cst{1,5}.visibleColor=[1 0 0];
cst{1,5}.visibleColor=[1 1 0];
cst{2,5}.visibleColor=[0 0 1];
cst{3,5}.visibleColor=[0 1 0];
cst{4,5}.visibleColor=[1 1 0];

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

if 0 % axis
idx=180:300;
idy=200:320;
idz=1:nz;
end

if 1 % sag
idx=150:360;
idy=1:ny;
idz=180:390;
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

matRadGUI
