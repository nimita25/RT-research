ptid='9306087'; 
%folder1 = ['G:\Tutorial_030624\code\' ptid '\'];
%folder2 = ['G:\Tutorial_030624\code\QC_' ptid '\'];
folder = ['..\' ptid '\'];
addpath('G:\Tutorial_030624\code\matRad-master\')
load([folder ptid '.mat'], 'ct', 'cst');
px = 2; % prescription dose
nfrac = 10; % number of fraction

ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc

ctv = cell(1,1);
ctv{1} = ctv1;
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);
    n_oar(i) = numel(oar{i});
end
c = [ctv; {body}; oar;];
%==============================================
N_c = numel(c);
n_c = zeros([N_c 1]);
for i = 1:N_c
    n_c(i) = numel(c{i});
end

name1=['Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-16-58.mat']; % Conv
name2=['Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-17-01.mat']; % GS
name3=['Results_9306087\res0202_9306087_4_NNZ_4_RND_50_2025-03-22-10-50.mat']; % 24 angles
name4=['Results_9306087\res0202_9306087_4_NNZ_4_RND_50_2025-03-22-10-52.mat']; % 72 angles


load(name4,'d');


d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,ct.cubeDim);
d(d>1.199)=1.199;
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([38],:);
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];%(!)


nx=ct.cubeDim(1);
ny=ct.cubeDim(2);
nz=ct.cubeDim(3);

idx=1:nx;
idy=1:ny;
idz=1:nz;



% for i=1:size(cst,1)
% mask=zeros([nx ny nz]);
% mask(cst{i,4}{1})=1;
% mask=mask(idx,idy,idz);
% cst{i,4}{1}=find(mask==1);
% end
matRadGUI