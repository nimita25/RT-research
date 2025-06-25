ptid='HN02'; 
folder = ['..\' ptid '\'];
addpath('G:\Tutorial_030624\code\utils')
addpath('G:\Tutorial_030624\code\matRad-master\')
load([folder ptid '.mat'], 'ct', 'cst');
px = 8; % prescription dose
nfrac = 5; % number of fraction

ctv1 = cst{15,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{9,4}{1};
N_oar = 4;
oar = cell(N_oar,1);
oar(1) = cst{2,4}(1); % R Parotid V50<30Gy
oar(2) = cst{11,4}(1); % OralCavity Dmean<40Gy
oar(3) = cst{17,4}(1); % Oropharynx Dmax<20Gy
oar(4) = cst{16,4}(1); % Larynx Dmax<20Gy

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

name1=['Results_1505_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-06-16-09-56.mat']; % Conv
name2=['Results_1505_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-06-16-09-57.mat']; % GS
name3=['Results_1505_HN02\res0202_HN02_72_NNZ_4_QC_50_2025-06-14-09-15.mat']; % 72 angles
name4=['Results_1505_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-06-19-11-29.mat']; % AG

% name1=['Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-11-37.mat']; % Conv
% name2=['Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-12-48.mat']; % 24 angles
% name3=['Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-11-39.mat']; % 72 angles

load(name4,'d');


d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,ct.cubeDim);
d(d>1.199)=1.199;
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([15],:);
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