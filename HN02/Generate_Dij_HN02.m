%clear;close all;clc

ptid='HN02'; % (!)
folder = ['D:\luoying\FLASH\re-irradiation cases\mat\HN02.mat'];
%addpath(genpath('D:\luoying\FLASH\matRad-master'));
%folder1 = ['D:\luoying\FLASH\IMPT\code' '\' ptid ];
% load([folder],'cst','ct');
load(['./' ptid '.mat'],'cst','ct');
addpath('./matRad-master/')

idctv=[15];
ctv0=cst{idctv(1),4}{1};

%bw = 5;  lss = 6;
s=3;r0=3;bw=5;lss=6; 

% ptv = ctv0;
% ptv0 = ptv;
ptv0=ctv2ptv_080720(ctv0,r0,ct.cubeDim,ct.resolution);
ptv=ctv2ptv_080720(ptv0,s,ct.cubeDim,ct.resolution);

[n1,n2] = size(cst);

cst_ptv = cell(n1+1,n2);
cst_ptv(1:n1,:) = cst;
cst_ptv(n1+1,:) = cst(idctv,:);
cst_ptv(n1+1,1) = {n1};
cst_ptv(n1+1,2) = {'PTV'};
cst_ptv(n1+1,3) = {'TARGET'};
cst_ptv{n1+1,4}{1} = ptv;
cst_ctv(n1+1,6) = {struct('type','square deviation','penalty',800,'dose',40,'EUD',NaN,'volume',NaN,'robustness','yes')};
for i = 1:n1
    cst_ptv(i,3) = {'OAR'};cst_ptv(i,6)={[]};
end

%% Plan Setup
pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth = bw;

%ID = [315,225];
ID = 0:15:345; %% delivery angle index

shift = [0,0,0];
% Dij0 = [];
% n = 0;
% id_angle = cell(numel(ID),1);
% Nray = cell(numel(ID),1);
% N_angle = numel(ID);
% N = zeros(numel(ID),1);
for i = 1:length(ID)
    id  = ID(i);
    angle = zeros(1,length(id));
    pln.propStf.gantryAngles    = id; % (!)
    pln.propStf.couchAngles     = angle; % (!)
    pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
    pln.propOpt.bioOptimization = 'none';
    pln.propOpt.runDAO          = false;
    pln.propOpt.runSequencing   = false;
    pln.propStf.isoCenter(1,:)
    
    stf = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
    stf.isoCenter = pln.propStf.isoCenter + shift;
    dij = matRad_calcParticleDose(ct,stf,pln,cst_ptv);
    %stf(i) = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
    %stf(i).isoCenter = pln.propStf.isoCenter + shift;
    %dij = matRad_calcParticleDose(ct,stf(i),pln,cst_ptv);
    %N(i) = size(dij.physicalDose{1},2);
    N = size(dij.physicalDose{1},2);
    %id_angle{i} = n + (1:N(i));
    %Nray{i} = stf(i).numOfBixelsPerRay;
    Nray = stf.numOfBixelsPerRay;
    %n = n + N(i);
    %Dij0 = [Dij0 dij.physicalDose{1}];
    save([ './' ptid '_' num2str(id)  '.mat'],'dij','stf', 'Nray','-v7.3');
end
%save([ folder1 '\' ptid '_v1'  '.mat'],'Dij0','stf','id_angle','N_angle','Nray','-v7.3');


% %% Set up uncertainty
% ss = 3;
% Shift = [[ss 0 0];[-ss 0 0];[0 ss 0];[0 -ss 0];[0 0 ss];[0 0 -ss]];
% for k = 1:6
%     Dij0 = [];
%     shift = Shift(k,:);
%     for i = 1:length(ID)
%         id  = ID(i);
%         angle = zeros(1,length(id));
%         pln.propStf.gantryAngles    = id; % (!)
%         pln.propStf.couchAngles     = angle; % (!)
%         % pln.propStf.gantryAngles    = [0 60 120 180 240]; % (!)
%         % pln.propStf.couchAngles     = angle; % (!)
%         pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
%         pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
%         pln.propOpt.bioOptimization = 'none';
%         pln.propOpt.runDAO          = false;
%         pln.propOpt.runSequencing   = false;
%         pln.propStf.isoCenter(1,:)
% 
%         stf(i) = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
%         stf(i).isoCenter = pln.propStf.isoCenter + shift;
%         dij = matRad_calcParticleDose(ct,stf(i),pln,cst_ptv);
%         Dij0 = [Dij0 dij.physicalDose{1}];
%     end
%     save([ folder1 '\' ptid '_v' num2str(k + 1)  '.mat'],'Dij0','stf','-v7.3');
% end
% 
% %% Range uncertainty
% c0 = 0.035;
% c = c0;
% shift = [0,0,0];
% Dij0 = [];
% for i = 1:length(ID)
%     id  = ID(i);
%     angle = zeros(1,length(id));
%     pln.propStf.gantryAngles    = id; % (!)
%     pln.propStf.couchAngles     = angle; % (!)
%     % pln.propStf.gantryAngles    = [0 60 120 180 240]; % (!)
%     % pln.propStf.couchAngles     = angle; % (!)
%     pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
%     pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
%     pln.propOpt.bioOptimization = 'none';
%     pln.propOpt.runDAO          = false;
%     pln.propOpt.runSequencing   = false;
%     pln.propStf.isoCenter(1,:)
% 
%     stf(i) = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
%     stf(i).isoCenter = pln.propStf.isoCenter + shift;
%     dij = matRad_calcParticleDose103118(1+c,ct,stf(i),pln,cst_ptv);
%     Dij0 = [Dij0 dij.physicalDose{1}];
% end
% save([ folder1 '\' ptid '_v' num2str(8)  '.mat'],'Dij0','stf','-v7.3');
% 
% c = -c0;
% shift = [0,0,0];
% Dij0 = [];
% for i = 1:length(ID)
%     id  = ID(i);
%     angle = zeros(1,length(id));
%     pln.propStf.gantryAngles    = id; % (!)
%     pln.propStf.couchAngles     = angle; % (!)
%     % pln.propStf.gantryAngles    = [0 60 120 180 240]; % (!)
%     % pln.propStf.couchAngles     = angle; % (!)
%     pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
%     pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
%     pln.propOpt.bioOptimization = 'none';
%     pln.propOpt.runDAO          = false;
%     pln.propOpt.runSequencing   = false;
%     pln.propStf.isoCenter(1,:)
% 
%     stf(i) = matRad_generateStf_112018_v2(lss,ct,cst_ptv,pln);
%     stf(i).isoCenter = pln.propStf.isoCenter + shift;
%     dij = matRad_calcParticleDose103118(1+c,ct,stf(i),pln,cst_ptv);
%     Dij0 = [Dij0 dij.physicalDose{1}];
% end
% save([ folder1 '\' ptid '_v' num2str(9)  '.mat'],'Dij0','stf','-v7.3');