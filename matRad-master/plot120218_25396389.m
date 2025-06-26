
clear

ptid='25396389';
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

if 1
load([folder 'Gmin5_mu0001_m70_n200_d.mat'],'d');
tmp=cst([3;2;1;17;15;],:);
cst=tmp;
else
load([folder 'robust_Gmin5_mu0001_m83_n200_v2_d.mat'],'d');
d=d{1};
tmp=cst([4;5;6;17;15;],:);
cst=tmp;
end

ctv1=cst{1,4}{1};   % ptv 53.9 (35*1.54)
ctv2=cst{2,4}{1};   % ptv 60.2 (35*1.72)
ctv3=cst{3,4}{1};   % ptv 70 (35*2)
ctv=cell(3,1);
ctv{1}=ctv1;
ctv{2}=ctv2;
ctv{3}=ctv3;
ctv_all=[];
for i=3:-1:1
    ctv{i}=setdiff(ctv{i},ctv_all);
    ctv_all=union(ctv_all,ctv{i});
end
cst{1,4}{1}=ctv{1};
cst{2,4}{1}=ctv{2};
cst{3,4}{1}=ctv{3};
cst{4,4}{1}=setdiff(cst{4,4}{1},ctv_all);
cst{5,4}{1}=setdiff(cst{5,4}{1},ctv_all);

d=interpn(x0,y0,z0,d,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d/2;


cst{3,5}.visibleColor=[0 0 0];
cst{2,5}.visibleColor=[1 0 1];
cst{1,5}.visibleColor=[1 1 1];
cst{4,5}.visibleColor=[0 0 1];
cst{5,5}.visibleColor=[0 1 0];

nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;

if 1 % axis
idx=130:270;
idy=145:285;
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
