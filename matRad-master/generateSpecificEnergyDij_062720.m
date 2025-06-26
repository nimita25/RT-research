
clear
close all
clc

ptid='51426287';

% addpath('C:\Users\hao\Desktop\Bowen\forBW_120519')
% Store dij in the specified folder
% ptid_1='C:\Users\hao\Desktop\adaptive grid\';

% if ~exist(ptid_1,'file')
%    mkdir(ptid_1);
% end

delta=5;dx=3;
% ptid='51426287'; % (!)
load(['C:\Users\hao\Desktop\FLASH\' ptid '.mat'],'cst','ct');
idctv=[9]; % (!)
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
% figure;imshow3D(mask,[]);
maskCube = mask;


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
% mask=zeros(ct.cubeDim);
% mask(ctv0)=mask(ctv0)+1;
% mask(ptv)=mask(ptv)+1;
% figure;imshow3D(mask,[]);
% 
% save([ptid '.mat'],'cst','ct','ptv');
pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
%% load machine
machine_fileName = [pln.radiationMode '_' pln.machine];
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false;
pln.propOpt.runSequencing   = false;
try
   load(machine_fileName);
   SAD = machine.meta.SAD;%? This subfield holds the geometrical source to axis distance [10000mm]
catch
   error(['Could not find the following machine file: ' machine_fileName ]); 
end
%% Specify a specific energy, here the maximum energy is used
Spec_energy = max([machine.data.energy]);

%% Specify specific bw, lss, sampling scale of angle and Number of angles
bw = 3; lss=3;  
N = 36;
% generate 36 gantryAngles
gantryAngles = 0 : 10 : (N-1)*10;
%%
[dij_cell,stf_cell] = generateMulti_dij(ptid,pln,ct,cst_ctv,cst_ptv,machine,Spec_energy,bw,lss,gantryAngles);

% Assemble multiple dij specified by sampling gantryAngles
% AssembleMode = 0;
% if AssembleMode
%    %% scale = 36,4,9,1;
%    % only angle 0
%    scale = 36;
%    [Dij,totalNumBixel] = aggregatedij(ptid_1,dij_cell,gantryAngles,scale,bw,lss);
%    
%    %including 9 angles
%    scale = 4;
%    [Dij,totalNumBixel] = aggregatedij(ptid_1,dij_cell,gantryAngles,scale,bw,lss);
%    
%    %including 4 angles
%    scale = 9;
%    [Dij,totalNumBixel] = aggregatedij(ptid_1,dij_cell,gantryAngles,scale,bw,lss);
%    
%    %including 36 angles
%    scale = 1;
%    [Dij,totalNumBixel] = aggregatedij(ptid_1,dij_cell,gantryAngles,scale,bw,lss);
% end


function [dij_cell,stf_cell] = generateMulti_dij(ptid,pln,ct,cst_ctv,cst_ptv,machine,Spec_energy,bw,lss,gantryAngles)
%%
% gantryAngles = 0 : 10 : (N-1)*10;
couchAngles  = zeros(1,numel(gantryAngles));
stf_cell = cell(1,numel(gantryAngles));
dij_cell = cell(1,numel(gantryAngles));
%%
for i = 1 : numel(gantryAngles)
    pln.propStfSpec_energy = Spec_energy;
    pln.numOfFractions  = 30;
    pln.propStf.bixelWidth      = bw;
    pln.propStf.longitudinalSpotSpacing = lss;

    pln.propStf.gantryAngles = gantryAngles(i);
    pln.propStf.couchAngles  = couchAngles(i); % (!)
    
    pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst_ptv,ct,0);
    pln.propStf.isoCenter(1,:)
    %%
%     [stf,~] = matRad_generateStf_maxenergy(ct,cst_ctv,pln,0);
    stf = matRad_generateStf_maxenergy270620(ct,cst_ctv,pln,0);
    stf_cell{i} = stf;
    %%
    shift=[0 0 0];
    stf2=stf;
    for i_1 = 1:pln.propStf.numOfBeams
        stf2(i_1).isoCenter=stf2(i_1).isoCenter+shift;
    end
    dij = matRad_calcParticleDose(ct,stf2,pln,cst_ctv);
    dij_cell{i} = dij;
    if ~isdeployed
%         save(['C:\Users\hao\Desktop\Bowen\51426287_33_080620\','51426287_dij_proton_',num2str(bw),num2str(lss),'_Spec_energy_',num2str(i),'.mat'],'dij','-v7.3')
%         save([ptid_1 '\','51426287_dij_proton_',num2str(bw),num2str(lss),'_Spec_energy_',num2str(i),'.mat'],'dij','-v7.3')
        save(['C:\Users\hao\Desktop\FLASH\' ptid '_dij_',num2str(i),'.mat'],'dij','-v7.3')
    end
end
%    %% Aggregated into a large Dij
%     if ~isdeployed
%         totalNumBixel = 0;
%         sampleAngle = downsample(gantryAngles,scale);
%         Id_sampleAgle = sampleAngle/10 + 1;
%         for j = Id_sampleAgle
% %             Dij = cat(2,Dij,dij_cell{j});
%             totalNumBixel = totalNumBixel + dij_cell{1, j}.totalNumOfBixels;
%             Dij = cat(2,Dij,dij_cell{1, j}.physicalDose{1});
%         end
%         save(['D:\图像处理\林博文\GH\New\forBW_120519\51426287\51426287_080620\','51426287_Dij_proton_',num2str(bw),num2str(lss),'_Spec_energy_NumAngle_',num2str(numel(gantryAngles)/scale),'.mat'],'Dij','-v7.3')
%     end

end


