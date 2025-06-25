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
% % Import rho values to calculate BED from physical dose
n_rho = numel(oar)+2;
rho = zeros([n_rho,1]);
rho(1) = cst{9,5}.betaX/cst{9,5}.alphaX;
rho(2) = cst{1,5}.betaX/cst{1,5}.alphaX;
rho(3) = cst{7,5}.betaX/cst{7,5}.alphaX;
rho(4) = cst{4,5}.betaX/cst{4,5}.alphaX;
rho(5) = cst{3,5}.betaX/cst{3,5}.alphaX;
rho_target = cst{9,5}.betaX/cst{9,5}.alphaX;
RBE = 1.1;
BED_target = (px*nfrac)*(RBE+rho_target*((px*nfrac)/nfrac));

K = 30; % number of iterations of the ADMM method
nm = [3 4 5];

% ==============================================================
%   N_obj     - Number of objectives.
%   w_obj     - Objective weight.
%   s_obj     - Prescription dose.
%   n_obj     - Array stores the number of active index for each objective.
%   c_obj     - Strucuture index for each objective.
%   id_obj    - Active index for each DVH objective.
%   type_obj  - DVH type of each objective.
%   0-> b = px/0 in the constr; 
%   1-> BED max constraint; 
%   2-> BED mean constraint; 
%   3-> BED DVH max constraint; 
%   4-> BED DVH min constraint for target; 
%   5-> Dmax constraint; 6-> D mean constraint; 
%   7-> DVH max constraint; 
%   8-> DVH min constraint
%   For target: type obj = [0;5;8]
% ==============================================================

N_obj = 10;
c_obj = [1; 1; 1; 3; 3; 3; 4; 5; 4; 5;];
type_obj = [0;5;8;2;3;0; 3; 2; 0; 0;];
s_obj = [px;1.1*px;px;18;12;0;27;20;0;0;];
n_obj = [nan;nan;0.95;nan;0.3;nan;0.6;nan;nan;nan;];
w_obj = [1;1;1;1;1;0.05;1;1;0.01;0.01;];

% N_obj = 12;
% type_obj = [0;2;3;0;3;1;1;0;1;1;0;0];
% w_obj = [1;1;1;0.1;1;1;1;0.005;1;1;0.01;0.01]; 
% s_obj = [px;px;px*1.1;0;px;[18;12;0;27;20;0;0]/nfrac];
% n_obj = round([nan;n_c(1)*0.95;nan;nan;nan;n_c(3)*0.5;n_c(3)*0.3;nan;n_c(4)*0.6;n_c(5)*0.5;nan;nan;]);
% c_obj = [1;1;1;2;2;3;3;3;4;5;4;5];


%% Set parameters for the Augmented Lagrangian -> check this!!!
mu_n0 = 1e-5;%0.1; %coefficient of penalty for linear constraint: dose to tumor <= 1.1px
mu_min = 1e-5;%0.1; %coefficient of penalty for DVH min constraint for tumor
mu = 1e-5;%0.1; %coefficient of penalty for BED max, BED mean constraint for OAR
mu_max = 1e-5;%0.1; %coefficient of penalty for BED DVH max constraint for OAR
mu_g = 1e-5;%0.01; %coefficient of penalty for the constraint: u >= 0
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

n_c = zeros([numel(c) 1]); %number of voxels in each DVH-max, DVH-mean OAR
for i = 1:numel(n_c)
    n_c(i) = numel(c{i});
end

% Update s_obj according to BED values
for i = 1:N_obj
    if type_obj(i) == 1 || type_obj(i) == 2 || type_obj(i) == 3
        s_obj(i) = s_obj(i)*(RBE+rho(c_obj(i))*(s_obj(i)/nfrac));
    end
end

out_folder = "output/";
resd = cell(3,1);
% fnames = {'./output/BED-ADMM-2024-09-02-08-22.mat' ...
%     './output/BED-ADMM-2024-09-02-10-09.mat',...
%     './output/BED-ADMM-2024-09-02-06-36.mat'};
fnames = {'./output/BED-ADMM-2024-09-02-06-36.mat' ...
    './output/BED-ADMM-2024-10-13-18-25.mat'};%'./output/BED-ADMM-2024-09-02-08-22.mat'};

% Load Dij matrices
unique_fields = 6;
% fname_Dij = strcat("Cost_matrix_",int2str(unique_fields),".mat");
% disp("Loading data from file...")
% load(fname_Dij);

% Initialize values for DVH plots
n_s = 100;
BED_T = zeros(n_s, numel(fnames));
BED_O = cell(numel(nm),1);
for nn = 1:numel(nm)
    BED_O{nn} = zeros(n_s, numel(fnames));
end

for ff = 1:numel(fnames)
    if ff == 1
        unique_fields = 1;
        % fname_Dij = strcat("Cost_matrix_",int2str(unique_fields),".mat");
        % disp("Loading data from file...")
        % load(fname_Dij);
    end
    % Load output here
    fname = fnames{ff};
    disp('Loading output...')
    load(fname);
    % x0 = max(0,fluencevector);
    % 
    % %% Normalize the dose
    % NF = zeros(min(nfrac,unique_fields),1);
    % n = 0;
    % d = zeros(size(Cost_matrix{1},1),min(nfrac,unique_fields));
    % for nf = 1:min(nfrac,unique_fields)
    %     [tmp_index1,tmp_index] = size(Cost_matrix{nf});
    %     y = Cost_matrix{nf}*x0(n+ (1:tmp_index));
    %     y = y(c{1});
    %     y2 = sort(y,'descend');
    %     NF(nf) = px/y2(ceil(n_c(1)*0.95));
    %     %d_tmp = Cost_matrix{nf}*x0(n+ (1:tmp_index))*NF(nf);
    %     %Dmax_tmp = max(Dmax_tmp,max(d_tmp(:))/px);
    %     d(:,nf) = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF(nf);
    %     n = n+tmp_index;
    % end
    

    % Get mean doses for further calculations
    TD = zeros(size(d,1),1);
    unique_fields = size(d,2);
    for nf = 1:nfrac
        TD = TD+d(:,mod(nf-1,unique_fields)+1);
    end
    %TD = sum(d,2);
    MD = TD/nfrac;
    BED2 = zeros(size(d,1),1);
    for nf = 1:nfrac
        BED2 = BED2+RBE*d(:,mod(nf-1,unique_fields)+1)+rho(3)*d(:,mod(nf-1,unique_fields)+1).^2;
    end
    
    % Get output metrics
    [D95, Dmax, CI, BEDmean_lung, BED30_lung, BEDmean_heart, BED60_heart, ...
    BEDmean_eso] =  calcparaBED_7119049(MD, BED2, px, c);
    resd{ff} = [D95, Dmax, CI, BEDmean_lung, BED30_lung, BEDmean_heart, BED60_heart, ...
    BEDmean_eso];

    % Generating DVH plot for target
    %DVH_eval = BED2(c{1})/BED_target;
    DVH_eval = MD(c{1})/px;
    N = numel(DVH_eval);
    % Create dose-volume histogram for target BED
    t = linspace(0, 1.2, n_s);
    for i = 1:n_s
        BED_T(i,ff) = numel(find(DVH_eval >= t(i))) / N;
    end

    % Generating DVH plot for OAR
    for nn = 1:numel(nm)
        DVH_eval = BED2(c{nm(nn)})/BED_target;
        N = numel(DVH_eval);
        % Create dose-volume histogram for target BED
        t = linspace(0, 1.2, n_s);
        for i = 1:n_s
            BED_O{nn}(i,ff) = numel(find(DVH_eval >= t(i))) / N;
        end
    end


end

name = strcat('output/metrics-',string(datetime('now','Format',"yyyy-MM-dd-HH-mm")),'.mat');
clear outp;
outp.resd = resd;
outp.fnames = fnames;
save(name,"-struct","outp");

% Generate plots
figure
hold on
t = linspace(0, 1.2, n_s);
tid = 100;
for ff = 1:numel(fnames)
    if ff == 1
    plot(t(1:tid)*100,100*BED_T((1:tid),ff),'LineWidth',2,'Color','red');
    else
        plot(t(1:tid)*100,100*BED_T((1:tid),ff),'LineWidth',2,'Color','blue');
    end
end
axis([90 110 0 100])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('ARC','multi-IMPT','Location','northeast','FontSize',14);
%legend('NC-BED','NC-Dose','Conv','Location','SouthWest');
grid on
%title(strcat('BED DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'BEDTumor_2908_',ptid,'.fig'));
hold off


figure
hold on
t = linspace(0, 1.2, n_s);

try
    % plot(t*100,BED_O{1}(:,1),'LineWidth',2,'LineStyle','--','Color','red','DisplayName','Lung');
    % plot(t*100,BED_O{1}(:,2),'LineWidth',2,'LineStyle','--','Color','blue','HandleVisibility','off');
    L(1:2) = plot(t*100,100*BED_O{1}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
    L(3:4) = plot(t*100,100*BED_O{1}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
catch
    
end
try
    % plot(t*100,BED_O{2}(:,1),'LineWidth',2,'LineStyle','-','Color','red','DisplayName','Heart');
    % plot(t*100,BED_O{2}(:,2),'LineWidth',2,'LineStyle','-','Color','blue','HandleVisibility','off');
    L(5:6) = plot(t*100,100*BED_O{2}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
    L(7:8) = plot(t*100,100*BED_O{2}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
catch
    
end
try
    % plot(t*100,BED_O{3}(:,1),'LineWidth',2,'LineStyle',':','Color','red','DisplayName','Esophagus');
    % plot(t*100,BED_O{3}(:,2),'LineWidth',2,'LineStyle',':','Color','blue','HandleVisibility','off');
    L(9:10) = plot(t*100,100*BED_O{3}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
    L(11:12) = plot(t*100,100*BED_O{3}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
catch
    
end
% for ff = 1:numel(fnames)
%     if ff == 1
%     plot(t,BED_O(:,ff),'LineWidth',2,'LineStyle',':','Color','red');
%     else
%         plot(t,BED_O(:,ff),'LineWidth',2,'Color','red');
%     end
%     %plot(t,BED_O(:,ff),'LineWidth',2);
% end
axis([0 110 0 80])
xlabel('BED (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
%legend('Lung','Heart','Esophagus','Location','northeast','FontSize',14);
% legend show
legend(L([2,6,10]), 'Lung','Heart','Esophagus','Location','northeast','FontSize',14)
grid on
%title(strcat('BED DVH in lung for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'BEDLung_2908_',ptid,'.fig'));
hold off



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Load output here
% disp('Loading output...')
% load('./output/BED-ADMM-2024-08-02-17-05.mat')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% objFn_track = zeros(K,1);
% Dmax_track = zeros(K,1);
% D95_track = zeros(K,1);
% 
% for i = 1:K
%     objFn_track(i) = violations{i}(1);
%     Dmax_track(i) = violations{i}(2);
%     D95_track(i) = violations{i}(3);
% end
% 
% x0 = max(0,fluencevector);
% d = zeros(size(Cost_matrix{1},1),1);
% 
% %% Normalize the dose
% NF = zeros(N_obj,1);
% n = 0;
% Dmax_tmp = 0;
% for nf = 1:N_obj
%     [tmp_index1,tmp_index] = size(Cost_matrix{mod(nf-1,unique_fields)+1});
%     y = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index));
%     y2 = sort(y,'descend');
%     NF(nf) = px/y2(ceil(tmp_index1*0.95));
%     d_tmp = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF(nf);
%     Dmax_tmp = max(Dmax_tmp,max(d_tmp(:))/px);
%     d = d+Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF(nf);
%     n = n+tmp_index;
% end
% d = d/N_obj;
% disp(Dmax_tmp);
% 
% %% Calculate Dmax, constraint values
% constr_values = cell(n_rho,1);
% for i = 1:n_rho
%     DVH_counter = 0;
%     constr_values{i} = zeros(length(type_constr{i}),1);
%     for ll = 1:length(type_constr{i})
%         if type_constr{i}(ll) == 1
%             tt = zeros(n_c(i),1);
%             n = 0;
%             for nf = 1:N_obj
%                 tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
%                 tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf);
%                 tt = tt+RBE*(tmp_prod)+rho(i)*(tmp_prod).^2;
%                 n = n+tmp_index;
%             end
%             constr_values{i}(ll) = max(tt,[],"all");
%         elseif type_constr{i}(ll) == 2
%             tt = 0;
%             n = 0;
%             for nf = 1:N_obj
%                 tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
%                 tt = tt+sum(RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf))+rho(i)*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf)).^2);
%                 n = n+tmp_index;
%             end
%             constr_values{i}(ll) = tt/n_c(i);
%         elseif type_constr{i}(ll) == 3
%             DVH_counter = DVH_counter+1;
%             DVH_eval = zeros(n_c(i),1);
%             n = 0;
%             for nf = 1:N_obj
%                 tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
%                 DVH_eval = DVH_eval + (RBE*tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf))+rho(i)*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf)).^2;
%                 n = n+tmp_index;
%             end
%             [eval2, I] = sort(DVH_eval,'descend');
%             ai = ceil(DVH_perc{i}(DVH_counter)*n_c(i));
%             constr_values{i}(ll) = eval2(ai);
%         end
%     end
% end
% 
% %% Calculate metrics for plot
% 
% 
% % Create dose-volume histogram
% d = d / px;
% N = numel(d);
% n = 100;
% disp(max(d(:)));
% 
% t = linspace(0, 1.1, n);
% dvh_O = zeros(n, 1);
% for i = 1:n
%     dvh_O(i) = numel(find(d >= t(i))) / N;
% end
% 
% DVH_eval = zeros(size(Cost_matrix{1},1),1);
% n = 0;
% for nf = 1:N_obj
%     tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
%     DVH_eval = DVH_eval + (RBE*Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF(nf))+rho_target*(Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF(nf)).^2;
%     n = n+tmp_index;
% end
% 
% DVH_eval = DVH_eval/BED_target;
% N = numel(DVH_eval);
% n = 100;
% 
% % Create dose-volume histogram for target BED
% t = linspace(0, 1.1, n);
% BED_O = zeros(n, 1);
% for i = 1:n
%     BED_O(i) = numel(find(DVH_eval >= t(i))) / N;
% end
% 
% % BED DVH for bladder
% nm = 1;  %1: bladder, 2: rectum, 3: femhead_lt, 4: femhead_rt, 5: penilebulb
% DVH_eval = zeros(size(tmp_Dij{1,nm},1),1);
% n = 0;
% for nf = 1:N_obj
%     tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
%     tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))*NF(nf);
%     DVH_eval = DVH_eval+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
%     n = n+tmp_index;
% end
% 
% DVH_eval = DVH_eval/BED_target;
% N = numel(DVH_eval);
% n = 100;
% 
% % Create dose-volume histogram for target BED
% t = linspace(0, 1.2, n);
% BED_BS_O = zeros(n, 1);
% for i = 1:n
%     BED_BS_O(i) = numel(find(DVH_eval >= t(i))) / N;
% end
% 
% %% Load data for conventional arc
% unique_fields = 1;
% for i = 1:N_obj %This can be changed as needed
%     %arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+18];
%     arc_id{i} = (0:23)';
% end
% if exist(strcat(out_folder,"Problem_para_conv_" ,ptid,".mat"),"file")
%     load(strcat(out_folder,"Problem_para_conv_" ,ptid,".mat"))
% else
%     Cost_matrix = cell([min(N_obj,unique_fields),1]);
%     tmp_Dij = cell(min(N_obj,unique_fields),n_rho);
% 
%     for i = 1:min(N_obj,unique_fields)
%         tmp_indices = [];
%         num_beamlets = zeros([numel(arc_id{i}),1]);
%         for j = 1:numel(arc_id{i})
%             disp('count');
%             disp(arc_id{i}(j));
%             load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
%             num_beamlets(j) = size(dij.physicalDose{1}, 2);
%         end
%         Cost_matrix{i} = sparse(numel(c{1}),sum(num_beamlets));
%         tmp_n = 0;
%         for m_c = 3:numel(c)
%             tmp_Dij{i,  m_c-2} = [];
%         end
%         for j = 1:numel(arc_id{i})
%             load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
%             disp('load');
%             disp(arc_id{i}(j));
%             Cost_matrix{i}(:,tmp_n + (1:num_beamlets(j))) = dij.physicalDose{1}(c{1},:); %Get rows corresponding to tumor voxels, and columns corresponding to arcs
%             tmp_n = tmp_n + num_beamlets(j);
%             for m_c = 3:numel(c)
%                 tmp_Dij{i,m_c-2} = [tmp_Dij{i,m_c-2},dij.physicalDose{1}(c{m_c},:)]; 
%             end
%         end
%     end
%     Cost_b = ones(numel(c{1}),1)*px;
%     dat.Cost_matrix = Cost_matrix;
%     dat.Cost_b = Cost_b;
%     dat.tmp_Dij = tmp_Dij;
%     fname = strcat(out_folder,'/Problem_para_conv_',ptid,'.mat');
%     if ~exist("output", 'dir')
%        mkdir("output")
%     end
%     save(fname,"-struct","dat");
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Load output here
% disp('Loading output...')
% load('./output/BED-ADMM-2024-06-30-21-45.mat')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% objFn_track_conv = zeros(K,1);
% Dmax_track_conv = zeros(K,1);
% D95_track_conv = zeros(K,1);
% 
% for i = 1:K
%     objFn_track_conv(i) = violations{i}(1);
%     Dmax_track_conv(i) = violations{i}(2);
%     D95_track_conv(i) = violations{i}(3);
% end
% 
% x0 = max(0,fluencevector);
% d = zeros(size(Cost_matrix{1},1),1);
% 
% n = 0;
% N_obj_c = 1;
% for nf = 1:N_obj_c
%     [tmp_index1,tmp_index] = size(Cost_matrix{mod(nf-1,unique_fields)+1});
%     y = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index));
%     y2 = sort(y,'descend');
%     NF = px/y2(ceil(tmp_index1*0.95));
%     d = d+Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF;
%     n = n+tmp_index;
% end
% 
% d = d / px;
% disp(max(d(:)));
% 
% %% Calculate Dmax, constraint values
% constr_values_conv = cell(n_rho,1);
% for i = 1:n_rho
%     DVH_counter = 0;
%     constr_values_conv{i} = zeros(length(type_constr{i}),1);
%     for ll = 1:length(type_constr{i})
%         if type_constr{i}(ll) == 1
%             tt = zeros(n_c(i),1);
%             n = 0;
%             for nf = 1:N_obj_c
%                 tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
%                 tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf);
%                 tt = tt+RBE*(tmp_prod)+rho(i)*(tmp_prod).^2;
%                 n = n+tmp_index;
%             end
%             tt = tt*N_obj;
%             constr_values_conv{i}(ll) = max(tt,[],"all");
%         elseif type_constr{i}(ll) == 2
%             tt = 0;
%             n = 0;
%             for nf = 1:N_obj_c
%                 tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
%                 tt = tt+sum(RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf))+rho(i)*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf)).^2);
%                 n = n+tmp_index;
%             end
%             tt = tt*N_obj;
%             constr_values_conv{i}(ll) = tt/n_c(i);
%         elseif type_constr{i}(ll) == 3
%             DVH_counter = DVH_counter+1;
%             DVH_eval = zeros(n_c(i),1);
%             n = 0;
%             for nf = 1:N_obj_c
%                 tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
%                 DVH_eval = DVH_eval + (RBE*tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf))+rho(i)*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))*NF(nf)).^2;
%                 n = n+tmp_index;
%             end
%             DVH_eval = DVH_eval*N_obj;
%             [eval2, I] = sort(DVH_eval,'descend');
%             ai = ceil(DVH_perc{i}(DVH_counter)*n_c(i));
%             constr_values_conv{i}(ll) = eval2(ai);
%         end
%     end
% end
% 
% % Create dose-volume histogram
% N = numel(d);
% n = 100;
% t = linspace(0, 1.1, n);
% dvh_C = zeros(n, 1);
% for i = 1:n
%     dvh_C(i) = numel(find(d >= t(i))) / N;
% end
% 
% DVH_eval = zeros(size(Cost_matrix{1},1),1);
% nf = 1;
% n = 0;
% tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
% DVH_eval = DVH_eval + (RBE*Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF)+rho_target*(Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF).^2;
% 
% DVH_eval = DVH_eval*nfrac;
% 
% DVH_eval = DVH_eval/BED_target;
% N = numel(DVH_eval);
% n = 100;
% 
% % Create dose-volume histogram for target BED
% t = linspace(0, 1.1, n);
% BED_C = zeros(n, 1);
% for i = 1:n
%     BED_C(i) = numel(find(DVH_eval >= t(i))) / N;
% end
% 
% % BED DVH for bladder
% %nm = 3; %1: bladder, 2: rectum, 3: femhead_lt, 4: femhead_rt, 5: penilebulb
% DVH_eval = zeros(size(tmp_Dij{1,nm},1),1);
% n = 0;
% for nf = 1:N_obj_c
%     tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
%     tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))*NF;
%     DVH_eval = DVH_eval+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
%     n = n+tmp_index;
% end
% 
% DVH_eval = DVH_eval*nfrac;
% 
% DVH_eval = DVH_eval/BED_target;
% N = numel(DVH_eval);
% n = 100;
% 
% % Create dose-volume histogram for target BED
% t = linspace(0, 1.2, n);
% BED_BS_C = zeros(n, 1);
% for i = 1:n
%     BED_BS_C(i) = numel(find(DVH_eval >= t(i))) / N;
% end
% 
% 
% %% Generate plots
% disp('Generating plots...')
% figure
% hold on
% plot(objFn_track);
% plot(objFn_track_conv);
% xlabel('K');
% ylabel('Objective function value');
% legend('Our method','Conventional');
% title(strcat('Objective function value vs K for case ID: ',ptid));
% saveas(gcf,strcat(out_folder,'ObjFnvK_C_',ptid,'.fig'));
% hold off
% 
% 
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
% 
% figure
% hold on
% n = 100;
% t = linspace(0, 1.1, n);
% plot(t, dvh_O, 'r', 'linewidth', 1);
% plot(t, dvh_C, 'b', 'linewidth', 1);
% xlabel('Dose (%)');
% ylabel('Volume (%)');
% legend('Our method','Conventional','Location','SouthWest');
% title(strcat('DVH in tumor for case ID: ',ptid));
% saveas(gcf,strcat(out_folder,'DVHTumor_C_',ptid,'.fig'));
% hold off
% 
% figure
% hold on
% plot(t, BED_O, 'r', 'linewidth', 1);
% plot(t, BED_C, 'b', 'linewidth', 1);
% xlabel('BED (%)');
% ylabel('Volume (%)');
% legend('Our method','Conventional','Location','SouthWest');
% title(strcat('BED DVH in tumor for case ID: ',ptid));
% saveas(gcf,strcat(out_folder,'BEDTumor_C_',ptid,'.fig'));
% hold off
% 
% 
% figure
% hold on
% n = 100;
% t = linspace(0, 1.2, n);
% plot(t, BED_BS_O, 'r', 'linewidth', 1);
% plot(t, BED_BS_C, 'b', 'linewidth', 1);
% xlabel('BED (%)');
% ylabel('Volume (%)');
% legend('Our method','Conventional');
% if nm == 5
%     title(strcat('BED DVH in penilebulb: ',ptid));
% elseif nm == 4
%     title(strcat('BED DVH in femhead-rt: ',ptid));
% elseif nm == 3 
%     title(strcat('BED DVH in femhead-lt: ',ptid));
% elseif nm == 2 
%     title(strcat('BED DVH in rectum: ',ptid));
% else
%     title(strcat('BED DVH in bladder: ',ptid));
% end
% %title(strcat('BED DVH in bladder: ',ptid));
% saveas(gcf,strcat(out_folder,'BEDBS_C_',ptid, '_nm', num2str(nm), '.fig'));
% hold off