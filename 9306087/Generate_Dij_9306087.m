% clear
% close all
% clc
% addpath(genpath('F:\KUMC\matRad-master'));
ptid = '9306087'; % (!)
folder = 'D:\Work\GH\LATTICE\Code_Energy\data\';
load(['./' ptid '.mat'],'cst','ct');
addpath('./matRad-master/')
% load([folder ptid '.mat'],'cst','ct');
idctv = [38]; % (!)
ctv0 = cst{idctv(1),4}{1};
ptv = ctv0;
ptv0 = ptv;
bw = 3;  lss = 6;
% ptv=ctv2ptv_080720(ctv0,r0,ct.cubeDim,ct.resolution);

[n1,n2] = size(cst);

cst_ptv = cell(n1+1,n2);
cst_ptv(1:n1,:) = cst;
cst_ptv(n1+1,:) = cst(idctv,:);
cst_ptv(n1+1,1) = {n1};
cst_ptv(n1+1,2) = {'PTV'};
cst_ptv(n1+1,3) = {'TARGET'};
cst_ptv{n1+1,4}{1} = ptv;
cst_ctv(n1+1,6) = {struct('type','square deviation','penalty',800,'dose',30,'EUD',NaN,'volume',NaN,'robustness','none')};
for i = 1:n1
    cst_ptv(i,3) = {'OAR'};cst_ptv(i,6)={[]};
end
% mask = zeros(ct.cubeDim);
% mask(ctv0) = mask(ctv0) + 1;
% mask(ptv) = mask(ptv) + 1;
% figure;imshow3D(mask,[]);

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 60;
pln.propStf.bixelWidth      = bw;

ID = 0:5:359; %% delivery angle index

shift = [0,0,0];
for i = 1:length(ID)
    id  = ID(i);
    disp(id)
    %angle = zeros(1,length(id));
    angle = 60;
    pln.propStf.gantryAngles    = id; % (!)
    pln.propStf.couchAngles     = angle; % (!)
    % pln.propStf.gantryAngles    = [0 60 120 180 240]; % (!)
    % pln.propStf.couchAngles     = angle; % (!)
    pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
    pln.propOpt.bioOptimization = 'none';
    pln.propOpt.runDAO          = false;
    pln.propOpt.runSequencing   = false;
    pln.propStf.isoCenter(1,:)
    
    stf = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
    stf.isoCenter = pln.propStf.isoCenter + shift;
    dij = matRad_calcParticleDose(ct,stf,pln,cst_ptv);
    save([ './' ptid '_' num2str(id) '_' num2str(angle)  '.mat'],'dij','stf','-v7.3');
end