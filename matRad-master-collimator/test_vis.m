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

% clear
% close all
% clc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
% load TG119.mat
% load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat
% load 14789903.mat % nose
% load 52721460.mat % nose
load 13440086.mat
% load 51584363.mat

ct.dicomInfo=[];


load resultGUI.mat

resultGUI=rmfield(resultGUI,{'overlapCube','physicalDose_beam1','physicalDose_beam2','w','wUnsequenced'});

load 13440086_dose.mat d_robust_proton d_robust_photon d_robust_hybrid
resultGUI.physicalDose=d_robust_proton{1};

tmp=cst([12;7;6;2],:);
% tmp=cst([12],:);
cst=tmp;

ct.cubeDim=[61 61 197];
idx=90:150;
idy=100:160;
ct.x=ct.x(idx);
ct.y=ct.y(idy);
ct.cubeHU{1}=ct.cubeHU{1}(idx,idy,:);
resultGUI.physicalDose=resultGUI.physicalDose(idx,idy,:);

for i=1:size(cst,1)
mask=zeros([260 260 197]);
mask(cst{i,4}{1})=1;
mask=mask(idx,idy,:);
cst{i,4}{1}=find(mask==1);
% cst{i,1}=i-1;
end

matRadGUI
