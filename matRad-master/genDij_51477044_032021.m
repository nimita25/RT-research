
clear
close all
clc
ptid='51477044'; % (!)

folder='C:\Users\hao\Desktop\matRad-master\';

load([folder ptid '\' ptid '.mat'],'cst','ct');
idctv=[10]; % (!)
ctv0=cst{idctv(1),4}{1};

[n1,n2]=size(cst);

% s=5;
s=0;
r=s+3;
bw=5;lss=6;

ctv=ctv2ptv_080720(ctv0,r,ct.cubeDim,ct.resolution);
cst_ctv=cell(n1+1,n2);
cst_ctv(1:n1,:)=cst;
cst_ctv(n1+1,:)=cst(idctv,:);
cst_ctv(n1+1,1)={n1};
cst_ctv(n1+1,2)={'CTV'};
cst_ctv{n1+1,4}{1}=ctv;
id=[n1+1];
for i=1:size(cst_ctv,1)
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

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.numOfFractions  = 30;
pln.propStf.bixelWidth      = bw;
pln.propStf.gantryAngles    = [0]; % (!)
pln.propStf.couchAngles     = [0]; % (!)
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

N=24;

gantryAngles=360/N*(0:(N-1));

for i=1:N
pln.propStf.gantryAngles = gantryAngles(i);
stf = matRad_generateStf_112018_v2(lss,ct,cst_ctv,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst);
save([folder ptid '\' ptid '_' num2str(gantryAngles(i)) '.mat'],'dij','stf','-v7.3');
end


return
shift=[0 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

shift=[-s 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

shift=[s 0 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

shift=[0 -s 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

shift=[0 s 0];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

shift=[0 0 -s];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

shift=[0 0 s];
stf2=stf;
for i=1:pln.propStf.numOfBeams
    stf2(i).isoCenter=stf2(i).isoCenter+shift;
end
dij = matRad_calcParticleDose(ct,stf2,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_' num2str(shift(1)) '_' num2str(shift(2)) '_' num2str(shift(3)) '.mat'],'dij','stf','-v7.3');

diff=-0.035;
dij = matRad_calcParticleDose103118(1+diff,ct,stf,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_range_m35.mat'],'dij','stf','-v7.3');
diff=+0.035;
dij = matRad_calcParticleDose103118(1+diff,ct,stf,pln,cst);
save([ptid '_r' num2str(r) '_bw' num2str(bw) '_lss' num2str(lss) '_range_35.mat'],'dij','stf','-v7.3');

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
