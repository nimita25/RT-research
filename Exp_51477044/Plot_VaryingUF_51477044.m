% ==============================================================
%%The lines of code here are unique for each test case
% ==============================================================
%% Define data path and parameters
addpath('../utils')
addpath('../utils')
ptid = '51477044';
folder = ['../' ptid '/'];
f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
t_px = 45; %total prescription dose
nfrac = 25; % number of fraction
px = t_px/nfrac; % prescription dose

%% Define target and OAR
ctv1 = cst{10,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{18,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{9,4}(1); % bladder D50: 25Gy; D20: 35Gy.
oar(2) = cst{17,4}(1); % rectum D50: 25Gy; D20: 35Gy; D10: 45Gy.
oar(3) = cst{11,4}(1); % femhead_lt D10: 25Gy.
oar(4) = cst{12,4}(1); % femhead_rt D10: 25Gy.
oar(5) = cst{14,4}(1); % penilebulb D50: 25Gy.
n_rho = numel(oar);
rho = zeros([n_rho,1]);
%rho(1) = cst{18,5}.betaX/cst{18,5}.alphaX;
rho(1) = cst{9,5}.betaX/cst{9,5}.alphaX;
rho(2) = cst{17,5}.betaX/cst{17,5}.alphaX;
rho(3) = cst{11,5}.betaX/cst{11,5}.alphaX;
rho(4) = cst{12,5}.betaX/cst{12,5}.alphaX;
rho(5) = cst{14,5}.betaX/cst{14,5}.alphaX;
rho_target = cst{10,5}.betaX/cst{10,5}.alphaX;
RBE = 1.1;

%% Define arc angles for each fraction
N_obj = nfrac;
arc_id = cell(N_obj,1); %each cell contains id of the arc associated with that objective
unique_fields = 6;
%unique_fields = 1;
for i = 1:N_obj %This can be changed as needed
    %arc_id{i} = [2*(i-1);2*(i+5)]; %;2*i-1;2*(i+5);2*(i+5)+1];
    %arc_id{i} = [6*(i-1)];
    %arc_id{i} = [2*(i-1);2*i-1;2*(i+5);2*(i+5)+1];
    arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+18];
    %arc_id{i} = (0:23)';
end
for i = 1:N_obj
    disp(arc_id{i}');
end

%% Define parameters of the constraints
%Type of constraint for each OAR
%type_constr = zeros([n_rho,1]);
type_constr = cell(n_rho,1); 
%type_constr = [3;3;3;3;3]; %1: BEDmax, 2: BEDmean, 3: BEDdvh
type_constr{1} = [3;3];
type_constr{2} = [3;3;3];
type_constr{3} = [3];
type_constr{4} = [3];
type_constr{5} = [3];

% DVH constraints for each OAR
num_DVH = [2;3;1;1;1];
DVH_perc = cell(n_rho,1);
DVH_perc{1} = [0.5;0.2];
DVH_perc{2} = [0.5;0.2;0.1];
DVH_perc{3} = 0.1;
DVH_perc{4} = 0.1;
DVH_perc{5} = 0.5;
DVH_max = cell(n_rho,1);
DVH_max{1} = [25;35];
DVH_max{2} = [25;35;45];
DVH_max{3} = 25;
DVH_max{4} = 25;
DVH_max{5} = 25;

%Dmax values of each OAR; used to calculate BEDmax; 0 -> no BEDmax constraint
Dmax = [0;0;0;0;0];

%Dmean values of each OAR; used to calculate BEDmean; 0 -> no BEDmean constraint
Dmean = [0;0;0;0;0];

%DVH min for target
DVH_min = px; %at least 95% of the voxels in tumor should get px dose in each fraction
DVH_min_perc = 0.95;

%% Set parameters for the Augmented Lagrangian
mu_n0 = 1e-5;%0.1; %coefficient of penalty for linear constraint: dose to tumor <= 1.1px
mu_min = 1e-5;%0.1; %coefficient of penalty for DVH min constraint for tumor
mu = 1e-5;%0.1; %coefficient of penalty for BED max, BED mean constraint for OAR
mu_max = 1e-5;%0.1; %coefficient of penalty for BED DVH max constraint for OAR
mu_g = 1e-5;%0.01; %coefficient of penalty for the constraint: u >= 0
K = 30; % number of iterations of the ADMM method
BED_target = (px*nfrac)*(RBE+rho_target*((px*nfrac)/nfrac));
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


UF_id = [12;8;6;4;3;1];
len_UF = length(UF_id);
objFn_track = zeros(K,len_UF);
Dmax_track = zeros(K,len_UF);
D95_track = zeros(K,len_UF);
n = 100;
%t = linspace(0, 1.1, n);
BED_O = zeros(n, len_UF);
BED_BS_O = zeros(n, len_UF);
dvh_O = zeros(n, len_UF);
tmp = [];


for uf = 1:length(UF_id)
    unique_fields = UF_id(uf);
    for i = 1:N_obj
        if unique_fields == 12
            arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+12];
        elseif unique_fields == 8
            arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+8;mod(i-1,unique_fields)+16];
        elseif unique_fields == 6
            arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+18];
        elseif unique_fields == 4
            arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+4;mod(i-1,unique_fields)+8;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+16;mod(i-1,unique_fields)+20];
        elseif unique_fields == 3
            arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+3;mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+9;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+15;mod(i-1,unique_fields)+18;mod(i-1,unique_fields)+21];
        elseif unique_fields == 2
            arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+2;mod(i-1,unique_fields)+4;mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+8;mod(i-1,unique_fields)+10;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+14;mod(i-1,unique_fields)+16;mod(i-1,unique_fields)+18;mod(i-1,unique_fields)+20;mod(i-1,unique_fields)+22];
        else
            arc_id{i} = (0:23)';
        end
    end
    disp("Loading data...")
    if exist(strcat(out_folder,"Problem_para_" ,ptid, "_UF", string(unique_fields), ".mat"),"file")
        load(strcat(out_folder,"Problem_para_" ,ptid, "_UF", string(unique_fields), ".mat"))
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
        fname = strcat(out_folder,"Problem_para_" ,ptid, "_UF", string(unique_fields), ".mat");
        %fname = strcat(out_folder,'/Problem_para_',ptid,'.mat');
        if ~exist("output", 'dir')
           mkdir("output")
        end
        save(fname,"-struct","dat");
    end
    if unique_fields == 12
        disp('Loading output...')
        load('./output/BED-ADMM-2024-07-30-13-08.mat')
    elseif unique_fields == 8
        disp('Loading output...')
        load('./output/BED-ADMM-2024-07-30-13-44.mat')
    elseif unique_fields == 6
        disp('Loading output...')
        load('./output/BED-ADMM-2024-07-30-14-49.mat')
    elseif unique_fields == 4
        disp('Loading output...')
        load('./output/BED-ADMM-2024-07-30-17-21.mat')
    elseif unique_fields == 3
        disp('Loading output...')
        load('./output/BED-ADMM-2024-07-31-12-09.mat')
    elseif unique_fields == 2
        disp('Loading output...')
        load('./output/BED-ADMM-2024-07-30-12-49.mat') %to update this later
    else
        disp('Loading output...')
        load('./output/BED-ADMM-2024-06-30-21-45.mat')
    end
    for i = 1:K
        objFn_track(i,uf) = violations{i}(1);
        Dmax_track(i,uf) = violations{i}(2);
        D95_track(i,uf) = violations{i}(3);
    end

    if unique_fields == 1
        N_obj = 1;
    end
    
    x0 = max(0,fluencevector);
    d = zeros(size(Cost_matrix{1},1),1);

    n = 0;
    for nf = 1:N_obj
        [tmp_index1,tmp_index] = size(Cost_matrix{mod(nf-1,unique_fields)+1});
        y = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index));
        y2 = sort(y,'descend');
        NF = px/y2(ceil(tmp_index1*0.95));
        d = d+Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF;
        n = n+tmp_index;
    end
    
    d = d/N_obj;
    d = d / px;
    N = numel(d);
    n = 100;
    
    % Create dose-volume histogram
    t = linspace(0, 1.1, n);
    for i = 1:n
        dvh_O(i) = numel(find(d >= t(i))) / N;
    end

    
    DVH_eval = zeros(size(Cost_matrix{1},1),1);
    n = 0;
    for nf = 1:N_obj
        tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
        DVH_eval = DVH_eval + (RBE*Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF)+rho_target*(Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF).^2;
        n = n+tmp_index;
    end
    if unique_fields == 1
        DVH_eval = DVH_eval*nfrac;
    end
    DVH_eval = DVH_eval/BED_target;
    N = numel(DVH_eval);
    n = 100;

    % Create dose-volume histogram for target BED
    t = linspace(0, 1.1, n);
    for i = 1:n
        BED_O(i,uf) = numel(find(DVH_eval >= t(i))) / N;
    end

    % BED DVH for bladder
    nm = 1; %1: bladder
    DVH_eval = zeros(size(tmp_Dij{1,nm},1),1);
    n = 0;
    for nf = 1:N_obj
        tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
        tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))*NF;
        DVH_eval = DVH_eval+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
        n = n+tmp_index;
    end

    if unique_fields == 1
        DVH_eval = DVH_eval*nfrac;
    end    
    %DVH_eval = DVH_eval/BED_max(nm);
    DVH_eval = DVH_eval/BED_target;
    N = numel(DVH_eval);
    n = 100;

    % Create dose-volume histogram for OAR BED
    t = linspace(0, 1.2, n);
    for i = 1:n
        BED_BS_O(i,uf) = numel(find(DVH_eval >= t(i))) / N;
    end

    % Legend name
    if unique_fields == 1
        tmp = [tmp; 'Conv'];
    else
        tmp = [tmp; string(24/unique_fields)+'fields/frac'];
    end


end



%% Generate plots
disp('Generating plots...')
figure
hold on
for uf = 1:len_UF
    plot(objFn_track(:,uf));
end
xlabel('K');
ylabel('Objective function value');
%legend('Our method','Conventional');
legend(tmp,'location','northeast','FontSize',10)
title(strcat('Objective function value vs K for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'ObjFnvK_C_',ptid,'.fig'));
hold off


% figure
% hold on
% plot(Dmax_track,'r');
% plot(D95_track,'g');
% plot(ones(K,1),'b');
% plot(Dmax_track_conv,'r','LineStyle','--');
% plot(D95_track_conv,'g','LineStyle','--');
% %plot(ones(K,1));
% xlabel('K');
% ylabel('Physical dose');
% legend('Dmax','D95','px');
% title(strcat('Dmax and D95 values in tumor vs K for case ID: ',ptid));
% saveas(gcf,strcat(out_folder,'DosevK_C_',ptid,'.fig'));
% hold off

% figure
% hold on
% n = 100;
% t = linspace(0, 1.1, n);
% for uf = 1:len_UF
%     plot(t, dvh_O(:,uf), 'linewidth', 1);
% end
% %plot(t, dvh_O, 'r', 'linewidth', 1);
% %plot(t, dvh_C, 'b', 'linewidth', 1);
% xlabel('Dose (%)');
% ylabel('Volume (%)');
% %legend('Our method','Conventional','Location','SouthWest');
% title(strcat('DVH in tumor for case ID: ',ptid));
% saveas(gcf,strcat(out_folder,'DVHTumor_C_',ptid,'.fig'));
% hold off

figure
hold on
for uf = 1:len_UF
    plot(t, BED_O(:,uf), 'linewidth', 1);
end
%plot(t, BED_O, 'r', 'linewidth', 1);
%plot(t, BED_C, 'b', 'linewidth', 1);
xlabel('BED (%)');
ylabel('Volume (%)');
legend(tmp,'location','southwest','FontSize',10)
%legend('Our method','Conventional','Location','SouthWest');
title(strcat('BED DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'BEDTumor_C_',ptid,'.fig'));
hold off


figure
hold on
n = 100;
t = linspace(0, 1.2, n);
for uf = 1:len_UF
    plot(t, BED_BS_O(:,uf), 'linewidth', 1);
end
%plot(t, BED_BS_O, 'r', 'linewidth', 1);
%plot(t, BED_BS_C, 'b', 'linewidth', 1);
xlabel('BED (%)');
ylabel('Volume (%)');
legend(tmp,'location','northeast','FontSize',10)
%legend('Our method','Conventional');
title(strcat('BED DVH in bladder: ',ptid));
saveas(gcf,strcat(out_folder,'BEDBS_C_',ptid,'.fig'));
hold off