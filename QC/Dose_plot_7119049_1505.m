ptid='7119049'; 
% folder1 = ['G:\Tutorial_030624\code\' ptid '\'];
% folder2 = ['G:\Tutorial_030624\code\QC_' ptid '\'];
folder = ['..\' ptid '\'];
addpath('G:\Tutorial_030624\code\matRad-master\')
load([folder ptid '.mat'], 'ct', 'cst');
px = 2; % prescription dose
nfrac = 30; % number of fraction

ctv1=cst{9,4}{1};   % ptv60, 2*30
body=cst{1,4}{1};
N_oar=3; % (!)
oar=cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
oar(4)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};

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

name1=['Results_1505_7119049\res0202_7119049_3_NNZ_3_RND_50_2025-06-16-11-49.mat']; % Conv
name2=['Results_1505_7119049\res0202_7119049_3_NNZ_3_RND_50_2025-06-16-12-29.mat']; % GS
name3=['Results_1505_7119049\res0202_7119049_72_NNZ_3_QC_50_2025-06-16-16-54.mat']; % 72 angles
name4=['Results_1505_7119049\res0202_7119049_3_NNZ_3_RND_50_2025-06-19-11-31.mat']; % AG

% name1=['Results_7119049\res0202_7119049_3_NNZ_3_RND_50_2025-03-22-12-17.mat']; % Conv
% name2=['Results_7119049\res0202_7119049_3_NNZ_3_RND_50_2025-03-22-13-41.mat']; % 24 angles
% name3=['Results_7119049\res0202_7119049_3_NNZ_3_RND_50_2025-03-22-13-55.mat']; % 72 angles

load(name1,'d');


d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,ct.cubeDim);
d(d>1.199)=1.199;
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([9],:);
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