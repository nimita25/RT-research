Tl = 35;
Td = 2;


ptid = '9306087';
folder = ['../' ptid '/'];
%f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc
n_rho = numel(oar)+2;
RBE = 1.1;
ctv = cell(1,1);
ctv{1} = ctv1; %row indices of Dij corresponding to target/tumor
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);  %row indices of Dij corresponding to OAR
    n_oar(i) = numel(oar{i}); %number of voxels in each OAR
end
c = [ctv; {body}; oar;]; %cell containing row indices corresponding to OAR, tumor, body

n_c = zeros([numel(c) 1]); %number of voxels in each DVH-max, DVH-mean OAR
for i = 1:numel(n_c)
    n_c(i) = numel(c{i});
end

% Load problem P1 data
% rhoT = 6, rhoOAR = 3
load('.\output\ST-2024-10-12-17-20.mat');
d = resT{12}.d;
nfrac = resT{12}.nfrac;
px = resT{12}.output_px;
y = d(c{1});
y2 = sort(y,'descend');
n_ctv95 = round(n_c(1)*0.95);
factor = px/y2(n_ctv95);
d = d*factor;

oarid = 3;
load('.\output\resADMM-2025-03-18-13-24.mat','d');
for nc = [oarid ]
    n_s = 100;
    DVHeval = d(c{nc})/px;
    N = numel(DVHeval);
    tt = linspace(0,1,n_s);
    v1 = zeros(n_s,1);
    for ii = 1:n_s
        v1(ii) = numel(find((DVHeval>=tt(ii))))/N;
    end
end

load('.\output\resADMM-2025-03-18-13-21.mat','d');
for nc = [oarid ]
    n_s = 100;
    DVHeval = d(c{nc})/px;
    N = numel(DVHeval);
    tt = linspace(0,1,n_s);
    v2 = zeros(n_s,1);
    for ii = 1:n_s
        v2(ii) = numel(find((DVHeval>=tt(ii))))/N;
    end
end

figure
hold on
plot(tt*100,v1*100,'LineWidth',2,'Color','red');
plot(tt*100,v2*100,'LineWidth',2,'Color','blue');
%axis([80 125 0 100])
legend('(P1)','(P1)+spatial opti','FontSize',14)
hold off