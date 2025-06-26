
clear
close all
clc
ptid='15581499'; % (!)
load([ptid '.mat'],'cst','ct');

body=cst{11,4}{1}(:);
ctv1=cst{23,4}{1}(:);
ctv2=cst{26,4}{1}(:);
ctv3=cst{13,4}{1}(:);
mask=zeros(ct.cubeDim);
mask(body)=mask(body)+1;
mask(ctv1)=mask(ctv1)+1;
mask(ctv2)=mask(ctv2)+1;
mask(ctv3)=mask(ctv3)+1;
figure;imshow3D(mask,[]);
ptv1=cst{24,4}{1}(:);
ptv2=cst{25,4}{1}(:);
ptv3=cst{22,4}{1}(:);
mask=zeros(ct.cubeDim);
mask(body)=mask(body)+1;
mask(ptv1)=mask(ptv1)+1;
mask(ptv2)=mask(ptv2)+1;
mask(ptv3)=mask(ptv3)+1;
figure;imshow3D(mask,[]);


id=[12;13;14;15;16;18;19;23;26;27;28;29;30;31;34;35;36;]; % (!)
cst(id,3)={'OAR'};cst(id,6)={[]};


pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [225 315 45 135]; % (!)
pln.propStf.couchAngles     = [0 0 0 0]; % (!)
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
save([ptid '_proton_0_0_0.mat'],'dij','-v7.3');
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
