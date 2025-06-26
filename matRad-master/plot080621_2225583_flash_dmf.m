
clear
pwd='C:\Users\hao\Desktop\FLASH_DMF_v2';
ptid='2225583';
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

% dx=3;
% shift=10;
% r=ceil(shift/dx);
% ctv1=cst{9,4}{1};
% 
% roi0=ctv2ptv(ctv1,r,ct.cubeDim);
% roi=setdiff(roi0,ctv1);
% cst{1,4}{1}=roi;    

px=15;
if 1
load([folder ptid '_px15_K4_mmu10_mup0.1_080421.mat'],'d','dr');
load([folder ptid '_px15_K4_mmu100_mup0.1_080421.mat'],'d','dr');
load([folder ptid '_px15_K4_mmu120_mup0.1_080421.mat'],'d','dr');
% load([folder ptid '_px15_K4_mmu140_mup0.1_mu0.001_mud1_080421.mat'],'d','dr');
if 1
d_c=8;
dr_c=40;
dmf=1.5;
load([folder ptid '_c.mat'],'c');
ctv=c{1};
var=struct('ctv',ctv,'dmf0',dmf,'d_c',d_c,'dr_c',dr_c,'k',100);
dd2=dr.*d;
alpha=dd2-dr_c*d;
beta=d-d_c;
[dmf,dmf_alpha,dmf_beta]=flash_dmf_logistic(alpha,beta,var);
de=d.*dmf;
d=de;
end

d=reshape(d,ctdim0)*1;
tmp=cst([53 55],:);

% tmp=cst([9 1],:);
% px=6;
cst=tmp;
cst{2,4}{1}=setdiff(cst{2,4}{1},cst{1,4}{1});

cst{1,5}.visibleColor=[1 1 0];
cst{2,5}.visibleColor=[0 0 1];
end

if 0
d_c=8;
dr_c=40;
dmf=1.5;

load([folder ptid '_c.mat'],'c');
ctv0=c{1};

% load([folder ptid '_px15_K3_mmu120_mup0.1_071621.mat'],'d','dr');
load([folder ptid '_px15_K3_mmu160_mup0.1_mu0.1_mud1_071621.mat'],'d','dr');

id=setdiff(find(dr>=dr_c&d>=d_c),ctv0);
de=d;
de(id)=d(id)/dmf;

d=reshape(de,ctdim0)*1;
tmp=cst([9],:);
% px=6;
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];
% cst{2,5}.visibleColor=[0 0 1];
end


if 0
d_c=8;
dr_c=40;
dmf=1.5;

load([folder ptid '_c.mat'],'c');
ctv0=c{1};

load([folder ptid '_px15_K3_mmu120_mup0.1_071621.mat'],'d','dr');
% load([folder ptid '_px15_K3_mmu160_mup0.1_mu0.1_mud1_071621.mat'],'d','dr');

id=zeros(ctdim0);
id(dr>=dr_c&d>=d_c)=1;
d=id;

tmp=cst([1],:);


px=1;
cst=tmp;
% cst{1,5}.visibleColor=[0 0 1];
cst{1,5}.visibleColor=[0 1 0];
% cst{2,5}.visibleColor=[0 1 0];
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
