
clear
close all
clc
ptid='51904269'; % (!)
folder='C:\Users\hao\Desktop\FLASH_v2\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
idctv=[36]; % (!)
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
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ctv,ct,0);
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false; 
pln.propOpt.runSequencing   = false;

N=36;gantryAngles=360/N*(0:(N-1));

for i=1:numel(gantryAngles)
pln.propStf.gantryAngles = gantryAngles(i);
stf = matRad_generateStf_112018_v2(lss,ct,cst_ctv,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst_ctv);
save([folder ptid '\' ptid '_' num2str(gantryAngles(i)) '_bp.mat'],'dij','stf','-v7.3');
end

machine_fileName = [pln.radiationMode '_' pln.machine];
load(machine_fileName);
Spec_energy = max([machine.data.energy]);

pln2=pln;
pln2.propStfSpec_energy = Spec_energy;

for i=1:numel(gantryAngles)
pln2.propStf.gantryAngles = gantryAngles(i);
stf = matRad_generateStf_tb(lss,ct,cst_ctv,pln2);
dij = matRad_calcParticleDose(ct,stf,pln,cst_ctv);
save([folder ptid '\' ptid '_' num2str(gantryAngles(i)) '_tb.mat'],'dij','stf','-v7.3');
end



