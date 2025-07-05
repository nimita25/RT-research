% clear
% close all
% clc

ptid='1243050'; 
gridDim = '111';
%load(['./' ptid '.mat'],'cst','ct');
addpath('../matRad-master-collimator/.')
addpath('../matRad-master/.')
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
%load([folder ptid '/' ptid '.mat'],'cst','ct');
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);
idctv=[11]; % (!) 

%% Generate ct, cst info according to grid dimension
[n1,n2]=size(cst);
for i=1:n1
tmpCube=zeros(ct.cubeDim);
tmpCube(cst{i,4}{1}) = 1;
% interpolate cube
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
cst{i,4}{1} = VdoseGrid;
end

x0=doseGrid.x;
y0=doseGrid.y;
z0=doseGrid.z;
x=ct.x;
y=ct.y;
z=ct.z;
[x0,y0,z0] = ndgrid(x0,y0,z0);
[x,y,z] = ndgrid(x,y,z);
ct.cubeHU{1}=interpn(x,y,z,ct.cubeHU{1},x0,y0,z0);
ct.resolution = doseGrid.resolution;
ct.x = doseGrid.x;
ct.y = doseGrid.y;
ct.z = doseGrid.z;
ct.cubeDim = doseGrid.cubeDim;


ctv0=cst{idctv(1),4}{1};

r0=3;
bw=5;lss=6;
ptv=ctv2ptv_080720(ctv0,r0,ct.cubeDim,ct.resolution);

a1=240;


cst_ptv=cell(n1+1,n2);
cst_ptv(1:n1,:)=cst;
cst_ptv(n1+1,:)=cst(idctv,:);
cst_ptv(n1+1,1)={n1};
cst_ptv(n1+1,2)={'PTV'};
cst_ptv(n1+1,3)={'TARGET'}; 
cst_ptv{n1+1,4}{1}=ptv;
cst_ctv(n1+1,6)={struct('type','square deviation','penalty',800,'dose',30,'EUD',NaN,'volume',NaN,'robustness','none')};
for i=1:n1
    cst_ptv(i,3)={'OAR'};cst_ptv(i,6)={[]};
end
mask=zeros(ct.cubeDim);
mask(ctv0)=mask(ctv0)+1;
mask(ptv)=mask(ptv)+1;
figure;imshow3D(mask,[]);

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = bw;
pln.propStf.gantryAngles    = [a1]; % (!)
pln.propStf.couchAngles     = [0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;
pln.propStf.isoCenter(1,:)

%%%%%%%%%%%% original dij generation
stf = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst_ptv);
%save([folder ptid '_proton_' num2str(a1) '.mat'],'dij','-v7.3');
save([folder '/dij_' ptid '_111_' num2str(pln.propStf.gantryAngles) '.mat'],'dij','stf','-v7.3');