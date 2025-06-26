
clear

ptid='15610140';
folder=[pwd '\' ptid '\'];


c_x=84;
c_y=84;
c_z=48;
nx0=256;
ny0=256;
nz0=32;

load([folder ptid '.mat'],'ct','cst');
x=ct.x;
y=ct.y;
z=ct.z;
ctdim=ct.cubeDim;
nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);


% load([folder ptid '_fine.mat'],'ct','cst');
% x=ct.x;
% y=ct.y;
% z=ct.z;

% ctdim=ct.cubeDim;

% [x0,y0,z0] = ndgrid(x0,y0,z0);
% [x,y,z] = ndgrid(x,y,z);

load([folder ptid '_d.mat'],'true_dose','pred_dose');
% dose1=dose;
% load([folder ptid '_ctvbodyoar_photon_setup_dose.mat'],'dose');
% dose2=dose;
% load([folder ptid '_res_photon_051820_d.mat'],'d2');
% dose2=reshape(d2,ctdim0);
% dose=0.7*dose1+0.3*dose2;
% dose=0.8*dose1+0.2*dose2;

d=zeros(ctdim);
% d(:,:,c_z+(1-nz0/2:nz0/2))=true_dose(nx0/2-c_x+(1:nx),ny0/2-c_y+(1:ny),:);
d(:,:,c_z+(1-nz0/2:nz0/2))=1/0.975*pred_dose(nx0/2-c_x+(1:nx),ny0/2-c_y+(1:ny),:);

cst=cst([16],:);
% d=dose;

px=1;



ctv=cst{1,4}{1};
% cst{2,4}{1}=setdiff(cst{2,4}{1},ctv);
% cst{3,4}{1}=setdiff(cst{3,4}{1},ctv);
% cst{4,4}{1}=setdiff(cst{4,4}{1},ctv);

% d=interpn(x0,y0,z0,d,x,y,z);
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d/px;


cst{1,5}.visibleColor=[1 1 0];
% cst{2,5}.visibleColor=[0 0 0];
% cst{3,5}.visibleColor=[0 0 1];
% cst{4,5}.visibleColor=[0 1 0];

% idx=1:nx;
% idy=1:ny;
% idz=1:nz;
% 
% if 0 % axis
% idx=60:300;
% idy=70:250;
% idz=1:nz;
% end
% 
% if 0 % sag
% idx=100:340;
% idy=1:ny;
% idz=99:339;
% end

% ct.cubeDim=[numel(idx) numel(idy) numel(idz)];
% 
% ct.x=ct.x(idx);
% ct.y=ct.y(idy);
% ct.z=ct.z(idz);
% ct.cubeHU{1}=ct.cubeHU{1}(idx,idy,idz);
% resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,idz);

% for i=1:size(cst,1)
% mask=zeros([nx ny nz]);
% mask(cst{i,4}{1})=1;
% mask=mask(idx,idy,idz);
% cst{i,4}{1}=find(mask==1);
% end

% figure;imagesc(resultGUI.physicalDose(:,:,1));colorbar;lim = caxis;caxis([0.0 1.2]*70);set(gca,'FontSize', 18);
matRadGUI
