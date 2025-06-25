% ==============================================================
%%The lines of code here are unique for each test case
% ==============================================================
%% Define data path and parameters
addpath('../utils')
ptid = '7119049';
folder = ['../' ptid '/'];
f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
t_px = 60; %total prescription dose
nfrac = 30; % number of fraction
px = t_px/nfrac; % prescription dose

%% Define target and OAR
ctv1 = cst{9,4}{1};   % ptv60, 2*30
body = cst{1,4}{1};
N_oar = 3;
oar = cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
n_rho = numel(oar);
rho = zeros([n_rho,1]);
%rho(1) = cst{18,5}.betaX/cst{18,5}.alphaX;
rho(1) = cst{7,5}.betaX/cst{7,5}.alphaX;
rho(2) = cst{4,5}.betaX/cst{4,5}.alphaX;
rho(3) = cst{3,5}.betaX/cst{3,5}.alphaX;
%rho(4) = cst{12,5}.betaX/cst{12,5}.alphaX;
%rho(5) = cst{14,5}.betaX/cst{14,5}.alphaX;
RBE = 1.1;

%% Define arc angles for each fraction
N_obj = nfrac;
arc_id = cell(N_obj,1); %each cell contains id of the arc associated with that objective
unique_fields = 6;
for i = 1:N_obj %This can be changed as needed
    %arc_id{i} = [2*(i-1);2*(i+5)]; %;2*i-1;2*(i+5);2*(i+5)+1];
    %arc_id{i} = [6*(i-1)];
    %arc_id{i} = [2*(i-1);2*i-1;2*(i+5);2*(i+5)+1];
    arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+18];
end
for i = 1:N_obj
    disp(arc_id{i}');
end

%% Define parameters of the constraints
%Type of constraint for each OAR
%type_constr = zeros([n_rho,1]);
type_constr = cell(n_rho,1); 
%type_constr = [3;3;3;3;3]; %1: BEDmax, 2: BEDmean, 3: BEDdvh
% type_constr{1} = [2;3];
% type_constr{2} = [3];
% type_constr{3} = [2];
type_constr{1} = [2;3];
type_constr{2} = [3];
type_constr{3} = [2];

% DVH constraints for each OAR
num_DVH = [1;1;0];
DVH_perc = cell(n_rho,1);
DVH_perc{1} = [0.3];
DVH_perc{2} = [0.6];
%DVH_perc{3} = 0.1;
DVH_max = cell(n_rho,1);
DVH_max{1} = [12];
DVH_max{2} = [27];
%DVH_max{3} = 25;

%Dmax values of each OAR; used to calculate BEDmax; 0 -> no BEDmax constraint
Dmax = [0;0;0];

%Dmean values of each OAR; used to calculate BEDmean; 0 -> no BEDmean constraint
Dmean = [18;0;20];

%DVH min for target
DVH_min = px; %at least 95% of the voxels in tumor should get px dose in each fraction
DVH_min_perc = 0.95;

%% Set parameters for the Augmented Lagrangian
mu_n0 = 1e-5;%0.1; %coefficient of penalty for linear constraint: dose to tumor <= 1.1px
mu_min = 1e-5;%0.1; %coefficient of penalty for DVH min constraint for tumor
mu = 1e-5;%0.1; %coefficient of penalty for BED max, BED mean constraint for OAR
mu_max = 1e-5;%0.1; %coefficient of penalty for BED DVH max constraint for OAR
mu_g = 1e-5;%0.01; %coefficient of penalty for the constraint: u >= 0
K = 50; % number of iterations of the ADMM method
% ==============================================================
%%End unique part of code for each test case
% ==============================================================

%%  Define optimization objective function parameters
% ==============================================================
%   c           - Row index for different structures.
%   Cost_matrix - Cost matrix for each fraction; includes beamlets of every
%                 field in the objective
%   Cost_b      - vector containing prescription dose for tumor voxel
% ==============================================================
ctv = cell(1,1);
ctv{1} = ctv1; %row indices of Dij corresponding to target/tumor
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);  %row indices of Dij corresponding to OAR
    n_oar(i) = numel(oar{i}); %number of voxels in each OAR
end
c = [ctv; {body}; oar;]; %cell containing row indices corresponding to OAR, tumor, body

n_c = zeros([n_rho 1]); %number of voxels in each DVH-max, DVH-mean OAR
for i = 1:n_rho
    n_c(i) = numel(c{i+2});
end

out_folder = "output/";

disp("Loading data...")
if exist(strcat(out_folder,"Problem_para_" ,ptid,".mat"),"file")
    load(strcat(out_folder,"Problem_para_" ,ptid,".mat"))
else

    Cost_matrix = cell([min(N_obj,unique_fields),1]);
    tmp_Dij = cell(min(N_obj,unique_fields),n_rho);
    
    for i = 1:min(N_obj,unique_fields)
        tmp_indices = [];
        num_beamlets = zeros([numel(arc_id{i}),1]);
        for j = 1:numel(arc_id{i})
            load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
            num_beamlets(j) = size(dij.physicalDose{1}, 2);
        end
        Cost_matrix{i} = sparse(numel(c{1}),sum(num_beamlets));
        tmp_n = 0;
        for m_c = 3:numel(c)
            tmp_Dij{i,  m_c-2} = [];
        end
        for j = 1:numel(arc_id{i})
            load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
            disp('load');
            disp(arc_id{i}(j));
            Cost_matrix{i}(:,tmp_n + (1:num_beamlets(j))) = dij.physicalDose{1}(c{1},:); %Get rows corresponding to tumor voxels, and columns corresponding to arcs
            tmp_n = tmp_n + num_beamlets(j);
            for m_c = 3:numel(c)
                tmp_Dij{i,m_c-2} = [tmp_Dij{i,m_c-2},dij.physicalDose{1}(c{m_c},:)]; 
            end
        end
    end
    Cost_b = ones(numel(c{1}),1)*px;
    dat.Cost_matrix = Cost_matrix;
    dat.Cost_b = Cost_b;
    dat.tmp_Dij = tmp_Dij;
    fname = strcat(out_folder,'/Problem_para_',ptid,'.mat');
    if ~exist("output", 'dir')
       mkdir("output")
    end
    save(fname,"-struct","dat");
end



%% Define parameters of the constraints
% ==============================================================
%   BED_max     - rhs of BED max constraints
%   BED_mean    - rhs of BED mean constraints
%   BED_DVH_max - rhs of BED DVH constraint
%   tmp_Dij     - coefficient matrix for all OAR
% ==============================================================
BED_max = zeros([n_rho,1]);
BED_mean = zeros([n_rho,1]);
BED_DVH_max = cell(n_rho,1);
for i = 1:n_rho
    for j = 1:length(type_constr{i})
        if type_constr{i}(j) == 1 
            BED_max(i) = Dmax(i)*(RBE+rho(i)*(Dmax(i)/nfrac));
        elseif type_constr{i}(j) == 2
            BED_mean(i) = Dmean(i)*(RBE+rho(i)*(Dmean(i)/nfrac));
        elseif type_constr{i}(j) == 3
            BED_DVH_max{i} = zeros(num_DVH(i),1);
            for n = 1:num_DVH(i)
                BED_DVH_max{i}(n) = DVH_max{i}(n)*(RBE+rho(i)*(DVH_max{i}(n)/nfrac));
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load output here
disp('Loading output...')
load('./output/BED-ADMM-2024-06-25-08-46.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

objFn_track = zeros(K,1);
Dmax_track = zeros(K,1);
D95_track = zeros(K,1);

for i = 1:K
    objFn_track(i) = violations{i}(1);
    Dmax_track(i) = violations{i}(2);
    D95_track(i) = violations{i}(3);
end


disp('Generating plots...')
figure
plot(objFn_track);
xlabel('K');
ylabel('Objective function value');
title(strcat('Objective function value vs K for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'ObjFnvK_',ptid,'.fig'));


figure
hold on
plot(Dmax_track);
plot(D95_track);
plot(ones(K,1));
xlabel('K');
ylabel('Physical dose');
legend('Dmax','D95','px');
title(strcat('Dmax and D95 values in tumor vs K for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DosevK_',ptid,'.fig'));
hold off

x0 = max(0,fluencevector);
d = zeros(size(Cost_matrix{1},1),1);

n = 0;
for nf = 1:N_obj
    tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
    d = d+Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index));
    n = n+tmp_index;
end

d = d/N_obj;

d = d / px;
N = numel(d);
n = 100;

% Create dose-volume histogram
t = linspace(0, 1.1, n);
dvh = zeros(n, 1);
for i = 1:n
    dvh(i) = numel(find(d >= t(i))) / N;
end
figure
plot(t, dvh, 'r', 'linewidth', 1);
xlabel('Dose (%)');
ylabel('Volume (%)');
title(strcat('DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DVHTumor_',ptid,'.fig'));




