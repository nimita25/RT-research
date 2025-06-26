% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
% load TG119.mat
% load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat
load 9984998.mat

cst{5,3}='OAR';cst{5,6}=[];
cst{10,3}='OAR';cst{10,6}=[];
cst{17,3}='OAR';cst{17,6}=[];
cst{18,3}='OAR';cst{18,6}=[];

% meta information for treatment plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon
% pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 30;

% beam geometry settings
% pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln.propStf.gantryAngles    = [0:72:359]; % [?]
% pln.propStf.couchAngles     = [0 0 0 0 0]; % [?]
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln.propStf.gantryAngles    = [210 260 310 0 50 100 150]; % [?]
% pln.propStf.couchAngles     = [0 0 0 0 0 0 0]; % [?]
pln.propStf.gantryAngles    = [90 270]; % [?]
pln.propStf.couchAngles     = [0 0]; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% shift=[1;1;]*[0 0 5];
% shift=[1;1;]*[100 100 100];
% pln.propStf.isoCenter=pln.propStf.isoCenter+shift;

% optimization settings
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% initial visualization and change objective function settings if desired
matRadGUI



%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
% stf = matRad_generateStf_101518(shift,ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    shift=[0 0 5];
    stf2=stf;
    for i=1:pln.propStf.numOfBeams
        stf2(i).isoCenter=stf2(i).isoCenter+shift;
    end
    dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
%     dij = matRad_calcParticleDose_101518(shift,ct,stf,pln,cst);
    shift=[0 0 0];
    stf2=stf;
    for i=1:pln.propStf.numOfBeams
        stf2(i).isoCenter=stf2(i).isoCenter+shift;
    end
%     dij = matRad_calcParticleDose(ct,stf2,pln,cst);
    r=+0.035;
    dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst);
end

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
% if strcmp(pln.radiationMode,'photons') && (pln.propertiesOpt.runSequencing || pln.propertiesOpt.runDAO)
if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

