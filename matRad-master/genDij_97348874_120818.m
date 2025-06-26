
clear
close all
clc
delta=5;dx=3;

ptid='97348874'; % (!)
load([ptid '.mat'],'cst','ct');
idctv=[19]; % (!)
ctv0=cst{idctv(1),4}{1};
[n1,n2]=size(cst);

ctv=ctv2ptv(ctv0,3,ct.cubeDim);
cst_ctv=cell(n1+1,n2);
cst_ctv(1:n1,:)=cst;
cst_ctv(n1+1,:)=cst(idctv,:);
cst_ctv(n1+1,1)={n1};
cst_ctv(n1+1,2)={'CTV'};
cst_ctv{n1+1,4}{1}=ctv;
id=[n1+1];
for i=1:size(cst,1)
    if ismember(i,id)
    cst_ctv(i,3)={'TARGET'};        
    else
    cst_ctv(i,3)={'OAR'};cst_ctv(i,6)={[]};
    end
end
mask=zeros(ct.cubeDim);
mask(ctv0)=mask(ctv0)+1;
mask(ctv)=mask(ctv)+1;
figure;imshow3D(mask,[]);

ptv=ctv2ptv(ctv0,ceil(delta/dx),ct.cubeDim);
cst_ptv=cell(n1+1,n2);
cst_ptv(1:n1,:)=cst;
cst_ptv(n1+1,:)=cst(idctv,:);
cst_ptv(n1+1,1)={n1};
cst_ptv(n1+1,2)={'PTV'};
cst_ptv{n1+1,4}{1}=ptv;
id=[n1+1];
for i=1:size(cst,1)
    if ismember(i,id)
    cst_ptv(i,3)={'TARGET'};        
    else
    cst_ptv(i,3)={'OAR'};cst_ptv(i,6)={[]};
    end
end
mask=zeros(ct.cubeDim);
mask(ctv0)=mask(ctv0)+1;
mask(ptv)=mask(ptv)+1;
figure;imshow3D(mask,[]);

save([ptid '.mat'],'cst','ct','ptv');

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [290 215 180.1]; % (!)
pln.propStf.couchAngles     = [0 0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

pln.propStf.isoCenter(1,:)

stf = matRad_generateStf_112018(ct,cst_ctv,pln);
% stf = matRad_generateStf(ct,cst_ctv,pln); %% DNU!!!!!!!!!!!!

shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_0_0.mat'],'dij','stf','-v7.3');
return

if 1
shift=[-delta 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_m' num2str(delta) '_0_0.mat'],'dij','-v7.3');
shift=[delta 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_' num2str(delta) '_0_0.mat'],'dij','-v7.3');
shift=[0 -delta 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_m' num2str(delta) '_0.mat'],'dij','-v7.3');
shift=[0 delta 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_' num2str(delta) '_0.mat'],'dij','-v7.3');
shift=[0 0 -delta];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_0_m' num2str(delta) '.mat'],'dij','-v7.3');
shift=[0 0 delta];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_0_' num2str(delta) '.mat'],'dij','-v7.3');
r=-0.05;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst_ctv);
save([ptid '_proton_range_m5.mat'],'dij','-v7.3');
r=+0.05;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst_ctv);
save([ptid '_proton_range_5.mat'],'dij','-v7.3');
end
% return

pln.radiationMode   = 'photons'; 
pln.propStf.gantryAngles    = [330 10 50 90 130 170 210]; % (!)
pln.propStf.couchAngles     = [0 0 0 0 0 0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
pln.propStf.isoCenter(1,:)

stf = matRad_generateStf(ct,cst_ptv,pln);

shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_0_0_0.mat'],'dij','-v7.3');
return
shift=[-delta 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_m' num2str(delta) '_0_0.mat'],'dij','-v7.3');
shift=[delta 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_' num2str(delta) '_0_0.mat'],'dij','-v7.3');
shift=[0 -delta 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_0_m' num2str(delta) '_0.mat'],'dij','-v7.3');
shift=[0 delta 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_0_' num2str(delta) '_0.mat'],'dij','-v7.3');
shift=[0 0 -delta];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_0_0_m' num2str(delta) '.mat'],'dij','-v7.3');
shift=[0 0 delta];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcPhotonDose(ct,stf2,pln,cst_ptv);
save([ptid '_photon_0_0_' num2str(delta) '.mat'],'dij','-v7.3');

