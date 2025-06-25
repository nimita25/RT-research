ptid='51477044'; 
folder = ['../' ptid '/'];
addpath('G:\Tutorial_030624\code\matRad-master\')
load([folder ptid '.mat'], 'ct', 'cst');
load('.\output\ST-2024-10-10-01-00.mat');
nfrac = resT{7}.nfrac;
px = resT{7}.output_px;


ctv1 = cst{10,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{18,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{9,4}(1); % bladder D50: 25Gy; D20: 35Gy.
oar(2) = cst{17,4}(1); % rectum D50: 25Gy; D20: 35Gy; D10: 45Gy.
oar(3) = cst{11,4}(1); % femhead_lt D10: 25Gy.
oar(4) = cst{12,4}(1); % femhead_rt D10: 25Gy.
oar(5) = cst{14,4}(1); % penilebulb D50: 25Gy.

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

name1=['.\output\resADMM-2025-03-18-13-30.mat']; % P1
name2=['.\output\resADMM-2025-03-18-13-31.mat']; % P1+spatial opti


load(name2,'d');


d=d/px;
%d(setdiff(1:prod(ct.cubeDim),c{2}))=0;
d=reshape(d,ct.cubeDim);
d(d>1.199)=1.199;
resultGUI=struct('physicalDose',[]);
resultGUI.physicalDose=d;

tmp=cst([10],:);
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