ptid='1243050'; 
folder = ['../' ptid '/'];
%addpath('G:\Tutorial_030624\code\matRad-master\')
load([folder ptid '.mat'], 'ct', 'cst');
px = 6; % prescription dose
nfrac = 4; % number of fraction

ctv1=cst{11,4}{1};   % ptv 24 Gy, 6 Gy*4
body=cst{1,4}{1};
N_oar=4;
oar=cell(N_oar,1);
oar(1)=cst{2,4}(1); % LargeBowel Dmax<38 Gy, V25<20 cc (D20cc<25 Gy)
oar(2)=cst{3,4}(1); % SmallBowel Dmax<35 Gy, V20<5 cc (D5cc<20 Gy)
oar(3)=cst{12,4}(1); % SpinalCord Dmax<25 Gy
oar(4)=cst{7,4}(1); % L_Kidney 150 cc<12 Gy (D150cc<12 Gy)
oar(5)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};

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

name1=['Results_1505_1243050/res0202_1243050_3_NNZ_3_RND_50_2025-06-16-09-46.mat']; % Conv
name2=['Results_1505_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-06-16-09-50.mat']; % GS
name3=['Results_1505_1243050\res0202_1243050_72_NNZ_3_QC_50_2025-06-14-16-29.mat']; % 72 angles
name4=['Results_1505_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-06-19-11-28.mat']; % AG
% name1=['Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-21-10-04.mat']; % Conv
% name2=['Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-21-09-53.mat']; % 24 angles
% name3=['Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-21-09-58.mat']; % 72 angles
% name1=['Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-02-24-10-41.mat']; % Conv
% name2=['Results_1243050\res0202_1243050_24_NNZ_3_QC_50_2025-03-03-09-26.mat']; % 24 angles
% name3=['Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-03-15-04.mat']; % 72 angles

load(name1,'d');


d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,ct.cubeDim);
d(d>1.199)=1.199;
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([11],:);
cst=tmp;
cst{1,5}.visibleColor=[1 1 0];%(!)
%disp(cst)

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