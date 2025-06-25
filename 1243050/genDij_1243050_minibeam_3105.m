clear;close all;clc;

ptid='1243050'; % (!)
%load(['C:\KUMC\mat\' ptid '.mat'],'cst','ct');
load(['./' ptid '.mat'],'cst','ct');
addpath('../matRad-master-collimator/.')
addpath('../matRad-master/.')
idctv=[11];
ctv0=cst{idctv(1),4}{1};

s=5;r0=3;bw=5;lss=6;

ptv0=ctv2ptv_080720(ctv0,r0,ct.cubeDim,ct.resolution);
ptv=ctv2ptv_080720(ptv0,s,ct.cubeDim,ct.resolution);

[n1,n2]=size(cst);

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
pln.propStf.gantryAngles = [240];
pln.propStf.couchAngles = [0];
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

stf = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);

%% add collimator
min_x = zeros(numel(stf), 1);
max_x = min_x;
min_y = min_x;
max_y = min_x;
for i = 1 : numel(stf)
    a = [stf(i).ray(:).rayPos_bev];
    min_x(i) = min(a(1:3:end));
    max_x(i) = max(a(1:3:end));
    min_y = min(a(3:3:end));
    max_y = max(a(3:3:end));
end
x_min = min(min_x);
x_max = max(max_x);
y_min = min(min_y);
y_max = max(max_y);

w_coll = [0.4, y_max-y_min];  % size of each open slit (mm)
cs_coll = 4;  % center-to-center distance (mm) 3-7
for ii = [0 2 1 -1 -0.5 0.5 1.5 -1.5]
cx = x_min : cs_coll : x_max;
shift = ii;
cx = cx+shift; %add shift to MSC position. Default is align MSC with x_min (left side)
ax = cx - w_coll(1);
bx = cx + w_coll(1);
ay = y_min * ones(size(ax));
by = y_max * ones(size(bx));

pln.collimator = [ax', bx', ay', by'];
pln.coll_iso = 200; % distance between collimator to isocenter (mm)
stf = getRatio_coll(stf, pln);

%% dose calculation
resolution=[1 1 3]; % mm
V = [cst{:,4}];
V = unique(vertcat(V{:}));
tmpCube    = zeros(ct.cubeDim);
tmpCube(V) = 1;
% interpolate cube
tmp_x=ct.resolution.x:resolution(1):ct.cubeDim(1)*ct.resolution.x;
tmp_y=(ct.resolution.y:resolution(2):ct.cubeDim(2)*ct.resolution.y)';
tmp_z=ct.resolution.z:resolution(3):ct.cubeDim(3)*ct.resolution.z;
DoseNumber = find(matRad_interp3((1:ct.cubeDim(1))*ct.resolution.x,(1:ct.cubeDim(2))*ct.resolution.y,(1:ct.cubeDim(3))*ct.resolution.z,...
    tmpCube,tmp_x,tmp_y,tmp_z,'nearest'));
cubeDim=[numel(tmp_x) numel(tmp_y) numel(tmp_z)];
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(cubeDim,DoseNumber);
xCoordsV = xCoordsV_vox(:)*resolution(1)-stf(1).isoCenter(1);
yCoordsV = yCoordsV_vox(:)*resolution(2)-stf(1).isoCenter(2);
zCoordsV = zCoordsV_vox(:)*resolution(3)-stf(1).isoCenter(3);
DosePoint  = [xCoordsV yCoordsV zCoordsV];
pln.DosePoint=DosePoint;

dij = matRad_calcParticleDose_point_coll(ct,stf,pln,cst_ptv);

Dij0=spalloc(prod(cubeDim),dij.totalNumOfBixels,1);
Dij0(DoseNumber,:)=dij.physicalDose{1};
dij.physicalDose{1}=Dij0;

save([ ptid '_' num2str(pln.propStf.gantryAngles) '_S' num2str(shift) '.mat'],'dij','-v7.3');
end