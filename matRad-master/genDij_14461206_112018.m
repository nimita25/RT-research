
clear
close all
clc
ptid='14461206'; % (!)
load([ptid '.mat'],'cst','ct');

id=[3;6;7;9;18]; % (!)
cst(id,3)={'OAR'};cst(id,6)={[]};

body=cst{21,4}{1}(:);
ctv1=cst{9,4}{1}(:);
ptv1=cst{14,4}{1}(:);

mask=zeros(ct.cubeDim);
mask(body)=mask(body)+1;
mask(ctv1)=mask(ctv1)+1;
mask(ptv1)=mask(ptv1)+1;

isempty(setdiff(ctv1,ptv1))

figure;imshow3D(mask,[]);
% figure;imshow3D(ct.cubeHU{1},[]);

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [195 270 345]; % (!)
pln.propStf.couchAngles     = [0 0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

stf = matRad_generateStf_112018(ct,cst,pln);

shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_0_0_0.mat'],'dij','stf','-v7.3');

% return

shift=[-5 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_m5_0_0.mat'],'dij','-v7.3');
shift=[5 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_5_0_0.mat'],'dij','-v7.3');
shift=[0 -5 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_0_m5_0.mat'],'dij','-v7.3');
shift=[0 5 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_0_5_0.mat'],'dij','-v7.3');
shift=[0 0 -5];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_0_0_m5.mat'],'dij','-v7.3');
shift=[0 0 5];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_0_0_5.mat'],'dij','-v7.3');
r=-0.035;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst);
save([ptid '_proton_range_m35.mat'],'dij','-v7.3');
r=+0.035;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst);
save([ptid '_proton_range_35.mat'],'dij','-v7.3');

return
pln.radiationMode   = 'photons'; 
pln.propStf.gantryAngles    = [210 260 310 0 50 100 150]; % (!)
pln.propStf.couchAngles     = [0 0 0 0 0 0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

stf = matRad_generateStf(ct,cst,pln);

shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_0_0_0.mat'],'dij','-v7.3');
shift=[-5 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_m5_0_0.mat'],'dij','-v7.3');
shift=[5 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_5_0_0.mat'],'dij','-v7.3');
shift=[0 -5 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_0_m5_0.mat'],'dij','-v7.3');
shift=[0 5 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_0_5_0.mat'],'dij','-v7.3');
shift=[0 0 -5];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_0_0_m5.mat'],'dij','-v7.3');
shift=[0 0 5];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
save([ptid '_photon_0_0_5.mat'],'dij','-v7.3');
