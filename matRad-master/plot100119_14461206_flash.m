
clear

ptid='14461206';
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

if 0
% load([folder '14461206_px2_K1_mmu140_NoOpt.mat'],'d');
% load([folder '14461206_px2_K1_mmu140.mat'],'d');
% load([folder '14461206_px6_K1_mmu140.mat'],'d');
% load([folder '14461206_px6_K3_mmu140.mat'],'d');
% load([folder '14461206_px6_K5_mmu125.mat'],'d');
% load([folder '14461206_px12_K1_mmu140.mat'],'d');
% load([folder '14461206_px12_K5_mmu130.mat'],'d');
% load([folder '14461206_px12_K9_mmu130.mat'],'d');
% load([folder '14461206_px20_K1_mmu140.mat'],'d');
load([folder '14461206_px20_K9_mmu120.mat'],'d');
load([folder '14461206_px20_K17_mmu120.mat'],'d');
d=reshape(d,ctdim0)*1;
tmp=cst([9],:);
px=20;
cst=tmp;
else
dx=3;
shift=5;
r=ceil(shift/dx);
ctv1=cst{9,4}{1};
roi=ctv2ptv(ctv1,r,ct.cubeDim);
cst{1,4}{1}=roi;    
    
    load([folder '14461206_px2_K1_mmu140_NoOpt.mat'],'dr');
load([folder '14461206_px2_K1_mmu140.mat'],'dr');
load([folder '14461206_px6_K1_mmu140.mat'],'dr');
load([folder '14461206_px6_K3_mmu140.mat'],'dr');
load([folder '14461206_px6_K5_mmu125.mat'],'dr');
load([folder '14461206_px12_K1_mmu140.mat'],'dr');
load([folder '14461206_px12_K5_mmu130.mat'],'dr');
load([folder '14461206_px12_K9_mmu130.mat'],'dr');
load([folder '14461206_px20_K1_mmu140.mat'],'dr');
load([folder '14461206_px20_K9_mmu120.mat'],'dr');
load([folder '14461206_px20_K17_mmu120.mat'],'dr');

d=reshape(dr,ctdim0);
tmp=cst([9 1],:);

px=1;
cst=tmp;
end


d=interpn(x0,y0,z0,d,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d/px;


% cst{2,4}{1}=setdiff(cst{2,4}{1},cst{1,4}{1});
% cst{3,4}{1}=setdiff(cst{3,4}{1},cst{1,4}{1});

cst{1,5}.visibleColor=[1 1 0];
cst{2,5}.visibleColor=[0 1 0];
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
