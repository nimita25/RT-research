
clear
close all
clc
ptid='52762512'; % (!)
load([ptid '.mat'],'cst','ct');

id=[1;27;28;29;30;32;37;42;44;45;46]; % (!)

cst(id,3)={'OAR'};cst(id,6)={[]};

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [10 350]; % (!)
pln.propStf.couchAngles     = [0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

stf = matRad_generateStf(ct,cst,pln);

shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_10_1.mat'],'dij','-v7.3');
shift=[-5*cos(10) -5*sin(10) 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_10_2.mat'],'dij','-v7.3');
shift=[5*cos(10) 5*sin(10) 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_10_3.mat'],'dij','-v7.3');
%shift=[0 -5 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
   % stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcParticleDose(ct,stf2,pln,cst);
%save([ptid '_proton_0_m5_0.mat'],'dij','-v7.3');
%shift=[0 5 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
   % stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcParticleDose(ct,stf2,pln,cst);
%save([ptid '_proton_0_5_0.mat'],'dij','-v7.3');
shift=[0 0 -5];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_10_4.mat'],'dij','-v7.3');
shift=[0 0 5];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_proton_10_5.mat'],'dij','-v7.3');
r=-0.035;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst);
save([ptid '_proton_10_6.mat'],'dij','-v7.3');
r=+0.035;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst);
save([ptid '_proton_10_7.mat'],'dij','-v7.3');

%pln.radiationMode   = 'photons'; 
%pln.propStf.gantryAngles    = [210 260 310 0 50 100 150]; % (!)
%pln.propStf.couchAngles     = [0 0 0 0 0 0 0]; % (!)
%pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
%pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%stf = matRad_generateStf(ct,cst,pln);

%shift=[0 0 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
    %stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_0_0_0.mat'],'dij','-v7.3');
%shift=[-5 0 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
    %stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_m5_0_0.mat'],'dij','-v7.3');
%shift=[5 0 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
    %stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_5_0_0.mat'],'dij','-v7.3');
%shift=[0 -5 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
   % stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_0_m5_0.mat'],'dij','-v7.3');
%shift=[0 5 0];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
    %stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_0_5_0.mat'],'dij','-v7.3');
%shift=[0 0 -5];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
   % stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_0_0_m5.mat'],'dij','-v7.3');
%shift=[0 0 5];
%stf2=stf;
%for i=1:pln.propStf.numOfBeams
   % stf2(i).isoCenter=stf2(i).isoCenter+shift;
%end
%dij = matRad_calcPhotonDose(ct,stf2,pln,cst);
%save([ptid '_photon_0_0_5.mat'],'dij','-v7.3');
