Tl = 35;
Td = 2;


ptid = '7119049';
folder = ['../' ptid '/'];
%f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
ctv1 = cst{9,4}{1};   % ptv60, 2*30
body = cst{1,4}{1};
N_oar = 3;
oar = cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
n_rho = numel(oar)+2;
RBE = 1.1;
% Calculate BED for target (rhs of equality constraint)
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
load('.\output\ST-2024-10-08-05-43.mat');
d = resT{7}.d;
nfrac = resT{7}.nfrac;
px = resT{7}.output_px;
y = d(c{1});
y2 = sort(y,'descend');
n_ctv95 = round(n_c(1)*0.95);
factor = px/y2(n_ctv95);
d = d*factor;

oarid = 1;
load('.\output\resADMM-2025-03-18-13-06.mat','d');
for nc = [oarid ]
    n_s = 100;
    DVHeval = d(c{nc})/px;
    N = numel(DVHeval);
    tt = linspace(0,1.2,n_s);
    v1 = zeros(n_s,1);
    for ii = 1:n_s
        v1(ii) = numel(find((DVHeval>=tt(ii))))/N;
    end
end

load('.\output\resADMM-2025-03-18-11-18.mat','d');
for nc = [oarid ]
    n_s = 100;
    DVHeval = d(c{nc})/px;
    N = numel(DVHeval);
    tt = linspace(0,1.2,n_s);
    v2 = zeros(n_s,1);
    for ii = 1:n_s
        v2(ii) = numel(find((DVHeval>=tt(ii))))/N;
    end
end

figure
hold on
plot(tt*100,v1*100,'LineWidth',2,'Color','red');
plot(tt*100,v2*100,'LineWidth',2,'Color','blue');
%axis([80 115 0 100])
%legend('spatially unoptimized','spatially optimized','FontSize',14)
hold off