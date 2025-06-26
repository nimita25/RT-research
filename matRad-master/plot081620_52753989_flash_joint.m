
clear
pwd='C:\Users\hao\Desktop\FLASH_joint';
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
% load([folder '\bw5\52753989_px6_mmu5_dro0_p98_mup0.1_mu0.1.mat'],'d');
load([folder '\bw5\52753989_px6_K3_mmu1240_mmu25_p98_s0.1_mup0.1_mu0.01.mat'],'d','d1','d2');
d=reshape(d2,ctdim0)*1;
tmp=cst([22],:);
px=6;
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];
else
dx=3;
shift=10;
r=ceil(shift/dx);
ctv1=cst{22,4}{1};

roi0=ctv2ptv(ctv1,r,ct.cubeDim);
roi1=setdiff(roi0,ctv1);
roi2=ctv2ptv(roi1,r,ct.cubeDim);
roi=intersect(roi0,roi2);

cst{1,4}{1}=roi;    
    
load([folder '\bw5\52753989_px6_mmu5_dro0_p98_mup0.1_mu0.1.mat'],'dr');
% load([folder '\bw5\52753989_px6_K3_mmu1240_mmu25_p98_s0.1_mup0.1_mu0.01.mat'],'dr');

d=reshape(dr,ctdim0);
tmp=cst([1],:);

px=1;
cst=tmp;
cst{1,5}.visibleColor=[0 1 0];
end


d=interpn(x0,y0,z0,d,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d/px;


% cst{2,4}{1}=setdiff(cst{2,4}{1},cst{1,4}{1});
% cst{3,4}{1}=setdiff(cst{3,4}{1},cst{1,4}{1});

% cst{2,5}.visibleColor=[0 0 1];
% cst{3,5}.visibleColor=[0 1 0];

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

% x=ones(nx,ny)*40;
% x(1)=100;
% figure;imagesc(resultGUI.physicalDose(:,:,1));colormap('jet');colorbar;lim = caxis;caxis([0 1.15]);set(gca,'FontSize', 18);
matRadGUI
