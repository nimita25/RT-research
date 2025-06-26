
clear
close all
clc
delta=5;dx=3;

ptid='9984998'; % (!)
load([ptid '.mat'],'cst','ct');
idctv=[10]; % (!)

[n1,n2]=size(cst);

cst_ctv=cst;
for i=1:n1
    if ismember(i,idctv)
    cst_ctv(i,3)={'TARGET'};        
    else
    cst_ctv(i,3)={'OAR'};cst_ctv(i,6)={[]};
    end
end

ctv=cst{idctv(1),4}{1};
ptv=ctv2ptv(ctv,ceil(delta/dx),ct.cubeDim);
cst_ptv=cell(n1+1,n2);
cst_ptv(1:n1,:)=cst;
cst_ptv(n1+1,:)=cst(idctv,:);
cst_ptv(n1+1,1)={n1};
cst_ptv(n1+1,2)={'PTV'};
cst_ptv{n1+1,4}{1}=ptv;
idptv=[n1+1];

for i=1:size(cst,1)
    if ismember(i,idptv)
    cst_ptv(i,3)={'TARGET'};        
    else
    cst_ptv(i,3)={'OAR'};cst_ptv(i,6)={[]};
    end
end

mask=zeros(ct.cubeDim);
mask(ctv)=mask(ctv)+1;
mask(ptv)=mask(ptv)+1;
figure;imshow3D(mask,[]);

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [90 270]; % (!)
pln.propStf.couchAngles     = [0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ctv,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

pln.propStf.isoCenter(1,:)

stf = matRad_generateStf_112018(ct,cst_ctv,pln);

shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_0_0.mat'],'dij','stf','-v7.3');

shift=[-delta 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_m5_0_0.mat'],'dij','-v7.3');
shift=[delta 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_5_0_0.mat'],'dij','-v7.3');
shift=[0 -delta 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_m5_0.mat'],'dij','-v7.3');
shift=[0 delta 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_5_0.mat'],'dij','-v7.3');
shift=[0 0 -delta];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_0_m5.mat'],'dij','-v7.3');
shift=[0 0 delta];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
save([ptid '_proton_0_0_5.mat'],'dij','-v7.3');
r=-0.035;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst_ctv);
save([ptid '_proton_range_m35.mat'],'dij','-v7.3');
r=+0.035;
dij = matRad_calcParticleDose103118(1+r,ct,stf,pln,cst_ctv);
save([ptid '_proton_range_35.mat'],'dij','-v7.3');


pln.radiationMode   = 'photons'; 
pln.propStf.gantryAngles    = [210 260 310 0 50 100 150]; % (!)
pln.propStf.couchAngles     = [0 0 0 0 0 0 0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ctv,ct,0);
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
save([ptid '.mat'],'cst','ct','ptv');

