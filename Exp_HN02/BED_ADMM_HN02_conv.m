% ==============================================================
%%The lines of code here are unique for each test case
% ==============================================================
%% Define data path and parameters
addpath('../utils')
ptid = 'HN02';
folder = ['../' ptid '/'];
f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
t_px = 40; %total prescription dose
nfrac = 5; % number of fraction
px = t_px/nfrac; % prescription dose

% ==============================================================
% Part of code unique for HN02
% ==============================================================
flag_dro = 1;
% ==============================================================
% Part of code unique for HN02
% ==============================================================

%% Define target and OAR
ctv1 = cst{15,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{9,4}{1};
N_oar = 4;
oar = cell(N_oar,1);
oar(1) = cst{2,4}(1); % R Parotid V50<30Gy
oar(2) = cst{11,4}(1); % OralCavity Dmean<40Gy
oar(3) = cst{17,4}(1); % Oropharynx Dmax<20Gy
oar(4) = cst{16,4}(1); % Larynx Dmax<20Gy
n_rho = numel(oar);
rho = zeros([n_rho,1]);
%rho(1) = cst{18,5}.betaX/cst{18,5}.alphaX;
rho(1) = cst{2,5}.betaX/cst{2,5}.alphaX;
rho(2) = cst{11,5}.betaX/cst{11,5}.alphaX;
rho(3) = cst{17,5}.betaX/cst{17,5}.alphaX;
rho(4) = cst{16,5}.betaX/cst{16,5}.alphaX;
rho_target = cst{15,5}.betaX/cst{15,5}.alphaX;
RBE = 1.1;

%% Define arc angles for each fraction
N_obj = nfrac;
arc_id = cell(N_obj,1); %each cell contains id of the arc associated with that objective
%unique_fields = 6;
unique_fields = 1;
for i = 1:N_obj %This can be changed as needed
    %arc_id{i} = [2*(i-1);2*(i+5)]; %;2*i-1;2*(i+5);2*(i+5)+1];
    %arc_id{i} = mod(i-1,unique_fields);%[6*(i-1)];
    %arc_id{i} = [2*(i-1);2*i-1;2*(i+5);2*(i+5)+1];
    %arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+18];
    arc_id{i} = (0:23)';
end
for i = 1:N_obj
    disp(arc_id{i}');
end

%% Define parameters of the constraints
%Type of constraint for each OAR
%type_constr = zeros([n_rho,1]);
wt_constr = cell(n_rho,1); 
wt_constr{1} = [1];
wt_constr{2} = [1];
wt_constr{3} = [1];
wt_constr{4} = [1];
%wt_constr{5} = [1];
type_constr = cell(n_rho,1); 
%type_constr = [3;3;3;3;3]; %1: BEDmax, 2: BEDmean, 3: BEDdvh
% type_constr{1} = [2;3];
% type_constr{2} = [3];
% type_constr{3} = [2];
type_constr{1} = [3];
type_constr{2} = [2];
type_constr{3} = [1];
type_constr{4} = [1];
%type_constr{5} = [3];


% DVH constraints for each OAR
num_DVH = [1;0;0;0];
DVH_perc = cell(n_rho,1);
%DVH_perc{1} = [0.3];
%DVH_perc{2} = [0.6];
%DVH_perc{3} = 0.1;
DVH_perc{1} = [0.5];  

DVH_max = cell(n_rho,1);
%DVH_max{1} = [12];
%DVH_max{2} = [27];
%DVH_max{3} = 25;
DVH_max{1} = 30;%0.12*t_px; %12

%Dmax values of each OAR; used to calculate BEDmax; 0 -> no BEDmax constraint
Dmax = [0;0;20;20];

%Dmean values of each OAR; used to calculate BEDmean; 0 -> no BEDmean constraint
Dmean = [0;40;0;0];

%DVH min for target
DVH_min = px; %at least 95% of the voxels in tumor should get px dose in each fraction
DVH_min_perc = 0.95;

%% Set parameters for the Augmented Lagrangian
mu_n0 = 1e-5;%0.1; %coefficient of penalty for linear constraint: dose to tumor <= 1.1px
mu_min = 1e-5;%0.1; %coefficient of penalty for DVH min constraint for tumor
mu = 1e-3;%0.1; %coefficient of penalty for BED max, BED mean constraint for OAR
mu_max = 1e-5;%0.1; %coefficient of penalty for BED DVH max constraint for OAR
mu_g = 1e-5;%0.01; %coefficient of penalty for the constraint: u >= 0
wt_obj = [1;1;40];
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

n_target = numel(c{1});
n_c = zeros([n_rho 1]); %number of voxels in each DVH-max, DVH-mean OAR
for i = 1:n_rho
    n_c(i) = numel(c{i+2});
end

Cost_matrix = cell([min(N_obj,unique_fields),1]);
tmp_Dij = cell(min(N_obj,unique_fields),n_rho);

for i = 1:min(N_obj,unique_fields)
    tmp_indices = [];
    num_beamlets = zeros([numel(arc_id{i}),1]);
    for j = 1:numel(arc_id{i})
        disp('count');
        disp(arc_id{i}(j));
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


% n = 0;
% tmp_Dij = cell(N_obj,n_rho);
% for m_c = 3:numel(c)
%     size_c = size(c{m_c},1);
%     for i = 1:N_obj
%         tmp_Dij{i,m_c-2} = [];
%         for j = 1:numel(arc_id{i})
%             load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
%             tmp_Dij{i,m_c-2} = [tmp_Dij{i,m_c-2},dij.physicalDose{1}(c{m_c},:)]; 
%         end
%     end
%     n=n+size_c;
% end


%% Define and initialize the decision variable 'u'. We use x0 in this code.
total_beamlets = 0;
N_obj_c = 1;
%for i = 1:N_obj
for i = 1:N_obj_c
    total_beamlets = total_beamlets+size(Cost_matrix{mod(i-1,unique_fields)+1},2);
    % for j = 1:numel(arc_id{i})
    %     load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
    %     total_beamlets = total_beamlets+size(dij.physicalDose{1}, 2);
    % end
end
%x0 = ones([total_beamlets,1])*10;
disp('initialize x');
x0 = zeros([total_beamlets,1]);


%% Initialize all variables in the ADMM method
y = x0;
%gamma = zeros(size(y));
gamma = y;
z_n0 = zeros(numel(c{1}),N_obj_c);
lambd_n0 = zeros(numel(c{1}),N_obj_c);
z_min = zeros(numel(c{1}),N_obj_c);
lambd_min = zeros(numel(c{1}),N_obj_c);
z_max = cell(n_rho,1);
lambd_max = cell(n_rho,1);
z_mean = cell(n_rho,1);
lambd_mean = cell(n_rho,1);
z_dvh = cell(n_rho,1);
lambd_dvh = cell(n_rho,1);
id = cell(n_rho,1); % set of active indices for each OAR
id_min = cell(N_obj_c,1); % set of active indices of DVH min in each fraction

n = 0;
disp('initialize z min');
for i = 1:N_obj_c
    tmp_index = size(Cost_matrix{mod(i-1,unique_fields)+1},2);
    z_n0(:,i) = Cost_matrix{mod(i-1,unique_fields)+1}*x0(n+ (1:tmp_index));
    z_min(:,i) = z_n0(:,i);
    n = n+tmp_index;
    [id_min{i},~] = f.find_active_indices_DVH_min(z_min(:,i),DVH_min_perc,Cost_matrix{mod(i-1,unique_fields)+1},DVH_min);
end

disp('initialize zmax');
for i = 1:n_rho
    for ll = 1:length(type_constr{i})
        if type_constr{i}(ll) == 1
            z_max{i} = zeros(n_c(i),N_obj_c);
            lambd_max{i} = zeros(n_c(i),N_obj_c);
            n = 0;
            for j = 1:N_obj_c
                tmp_index = size(tmp_Dij{mod(j-1,unique_fields)+1,i},2);
                z_max{i}(:,j) = tmp_Dij{mod(j-1,unique_fields)+1,i}*x0(n+ (1:tmp_index));
                %lambd_max{i}(:,j) = z{i}(:,j);
                n = n+tmp_index;
            end
        elseif type_constr{i}(ll) == 2
            z_mean{i} = zeros(n_c(i),N_obj_c);
            lambd_mean{i} = zeros(n_c(i),N_obj_c);
            n = 0;
            for j = 1:N_obj_c
                tmp_index = size(tmp_Dij{mod(j-1,unique_fields)+1,i},2);
                z_mean{i}(:,j) = tmp_Dij{mod(j-1,unique_fields)+1,i}*x0(n+ (1:tmp_index));
                %lambd_mean{i}(:,j) = z{i}(:,j);
                n = n+tmp_index;
            end
        elseif type_constr{i}(ll) == 3
            z_dvh{i} = zeros(n_c(i),N_obj_c);
            lambd_dvh{i} = zeros(n_c(i),N_obj_c);
            n = 0;
            for j = 1:N_obj_c
                tmp_index = size(tmp_Dij{mod(j-1,unique_fields)+1,i},2);
                z_dvh{i}(:,j) = tmp_Dij{mod(j-1,unique_fields)+1,i}*x0(n+ (1:tmp_index));
                %lambd_dvh{i}(:,j) = z{i}(:,j);
                n = n+tmp_index;
            end
        end
    end
end
disp('update ac');
for i = 1:n_rho
    %Update active index set for DVH max
    if any(type_constr{i} == 3)
        id{i} = cell(num_DVH(i),1);
        DVH_eval = zeros(n_c(i),1);
        n = 0;
        for nf = 1:N_obj_c
            tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
            DVH_eval = DVH_eval + (RBE*tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index)))+rho(i)*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))).^2;
            n = n+tmp_index;
        end
        DVH_eval = DVH_eval*N_obj;
        if strcmp(ptid,'9306087')
            [id{i},~] = f.find_active_indices_DVH_max(DVH_eval,num_DVH(i),DVH_perc{i}/n_c(i),n_c(i),BED_DVH_max{i},id{i});
        else
            [id{i},~] = f.find_active_indices_DVH_max(DVH_eval,num_DVH(i),DVH_perc{i},n_c(i),BED_DVH_max{i},id{i});
        end
    end
end


%% Define parameters to be saved after the experiment
outp.K = K;
outp.mu_names = ["mu_n0","mu_min","mu","mu_max","mu_g"];
outp.mu = [mu_n0,mu_min,mu,mu_max,mu_g];
outp.nfrac = nfrac;
outp.t_px = t_px;
outp.px = px;
outp.arc_id = arc_id;
outp.viol_names = ["Objective_value","Dmax","DVHmin","Viol_BEDmax","Viol_BEDmean","Viol_BEDDVH"];
outp.violations = cell(K,1);

%% Start the iterations of the ADMM method
for k = 1:K

    %% Update y
    %y = x0 - gamma;
    y = max(0,x0 - gamma);

    %% Update z_n0
    n = 0;
    disp('update z no');
    for nf = 1:N_obj_c
        tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
        z_n0(:,nf) = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index)) - lambd_n0(:,nf);
        n = n+tmp_index;
    end
    %ind = find(z_n0>1.1*px);
    %z_n0(ind) = 1.1*px;
    z_n0 = min(z_n0,1.1*px);

    %% Update z_min
    n = 0;
    disp('update z min');
    for nf = 1:N_obj_c
        tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
        if ~isempty(id_min{nf})
            z_min(id_min{nf},nf) = Cost_matrix{mod(nf-1,unique_fields)+1}(id_min{nf},:)*x0(n+ (1:tmp_index))- lambd_min(id_min{nf},nf);
            z_min(id_min{nf},nf) = max(z_min(id_min{nf},nf),DVH_min);
            % for i = 1:numel(id_min{nf})
            %     %z_min(id_min{nf}(i),nf) = Cost_matrix{mod(nf-1,unique_fields)+1}(id_min{nf}(i),:)*x0(n+ (1:tmp_index))- lambd_min(id_min{nf}(i),nf);
            %     if z_min(id_min{nf}(i),nf) < DVH_min
            %         z_min(id_min{nf}(i),nf) = DVH_min;
            %     end
            % end
        end
        n = n+tmp_index;
    end

    %% Update z^{nm}
    for nm = 1:n_rho
        for ll = 1:length(type_constr{nm})
            if type_constr{nm}(ll) == 1
                tmp_prod = zeros(n_c(nm),N_obj_c);
                %tt = zeros(n_c(nm,1));
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                    tmp_prod(:,nf) = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))-lambd_max{nm}(:,nf);
                    %tt = tt+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
                    %z_max{nm}(:,nf) = tmp_prod;
                    n = n+tmp_index;
                end
                %p = zeros(n_c(nm),1);
                %disp(size(tmp_prod));
                %disp(size(sum((tmp_prod+(RBE/(2*rho(nm)))).^2,2)));
                p = 1 - sqrt(((BED_max(nm)/(rho(nm)*N_obj))+(N_obj_c*RBE^2/(2*rho(nm))^2))./sum((tmp_prod+(RBE/(2*rho(nm)))).^2,2));
                % disp(min(p));
                % disp(max(p));
                p = max(0,p);
                %disp(size(p));
                z_max{nm} = (1-p).*tmp_prod-p.*RBE/(2*rho(nm));
                tmp = sum(z_max{nm},2)*RBE+rho(nm)*sum(z_max{nm}.^2,2);
                % disp(max(tmp));
                % disp(BED_max(nm));
                % tmp_id = find(tt>BED_max(nm));
                % for ii = 1:length(tmp_id)
                %     i = tmp_id(i);
                %     disp('optimizing Dmax');
                %     % Solve the optimization problem
                %     A = [];
                %     b = [];
                %     Aeq = [];
                %     beq = [];
                %     ceq = [];
                %     ub = [];
                %     lb = [];
                %     x = ones(N_obj,1)*1;
                %     obj_handle = @(x)f.define_BED_max_objective(x,N_obj,tmp_Dij,i,nm,lambd_max,x0,unique_fields);
                %     cineq_handle = @(x)f.define_BED_max_constraint(x,N_obj,nm,rho,BED_max(nm),RBE);
                %     options = optimoptions('fmincon','Display','off');
                %     [x,fval,exitflag,output] = fmincon(obj_handle,x,A,b,Aeq,beq,lb,ub,cineq_handle,options);
                %     z_max{nm}(i,:) = x;
                % end

                % for i = 1:n_c(nm)
                %     tt = 0;
                %     n = 0;
                %     for nf = 1:N_obj
                %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %         tt = tt+RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_max{nm}(i,nf))+rho(nm)*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_max{nm}(i,nf))^2;
                %         n = n+tmp_index;
                %     end
                %     if tt <= BED_max(nm)
                %         n = 0;
                %         for nf = 1:N_obj
                %             tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %             z_max{nm}(i,nf) = tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_max{nm}(i,nf);
                %             n = n+tmp_index;
                %         end
                %     else
                %         disp('optimizing Dmax');
                %         % Solve the optimization problem
                %         A = [];
                %         b = [];
                %         Aeq = [];
                %         beq = [];
                %         ceq = [];
                %         ub = [];
                %         lb = [];
                %         x = ones(N_obj,1)*1;
                %         obj_handle = @(x)f.define_BED_max_objective(x,N_obj,tmp_Dij,i,nm,lambd_max,x0,unique_fields);
                %         cineq_handle = @(x)f.define_BED_max_constraint(x,N_obj,nm,rho,BED_max(nm),RBE);
                %         options = optimoptions('fmincon','Display','off');
                %         [x,fval,exitflag,output] = fmincon(obj_handle,x,A,b,Aeq,beq,lb,ub,cineq_handle,options);
                %         z_max{nm}(i,:) = x;
                %     end
                % end
            elseif type_constr{nm}(ll) == 2
                %Update mean constraints
                disp('update Dmean')
                tmp_prod = zeros(n_c(nm),N_obj_c);
                %tt = zeros(n_c(nm,1));
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                    tmp_prod(:,nf) = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))-lambd_mean{nm}(:,nf);
                    %tt = tt+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
                    %z_max{nm}(:,nf) = tmp_prod;
                    n = n+tmp_index;
                end
                p = 1 - sqrt(((BED_mean(nm)*n_c(nm)/(N_obj*rho(nm)))+(n_c(nm)*N_obj_c*RBE^2/(2*rho(nm))^2))/sum((tmp_prod+(RBE/(2*rho(nm)))).^2,"all"));
                % disp(min(p));
                % disp(max(p));
                p = max(0,p);
                z_mean{nm} = ((1-p)*tmp_prod)-(p*RBE/(2*rho(nm)));
                tmp = (sum(z_mean{nm},"all")*RBE)+(rho(nm)*sum(z_mean{nm}.^2,"all"));
                % disp(max(tmp));
                % disp(n_c(nm)*BED_mean(nm));
                % tt = 0;
                % n = 0;
                % for nf = 1:N_obj
                %     tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %     tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))-lambd_mean{nm}(:,nf);
                %     z_mean{nm}(:,nf) = tmp_prod;
                %     tt = tt+sum(RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2);
                %     n = n+tmp_index;
                % end
                % for i = 1:n_c(nm)
                %     n = 0;
                %     for nf = 1:N_obj
                %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %         tt = tt+RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_mean{nm}(i,nf))+rho(nm)*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_mean{nm}(i,nf))^2;
                %         n = n+tmp_index;
                %     end
                % end
                % if tt <= n_c(nm)*BED_mean(nm)
                % 
                %     n = 0;
                %     for nf = 1:N_obj
                %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %         z_mean{nm}(:,nf) = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))-lambd_mean{nm}(:,nf);
                %         n = n+tmp_index;
                %     end
                % else
                %     disp('optimizing for BED mean')
                %     % Solve the optimization problem
                %     A = [];
                %     b = [];
                %     Aeq = [];
                %     beq = [];
                %     ceq = [];
                %     ub = [];
                %     lb = [];
                %     x = zeros(n_c(nm),N_obj);
                %     obj_handle = @(x)f.define_BED_mean_objective(x,N_obj,tmp_Dij,nm,lambd_mean,x0,n_c(nm),unique_fields);
                %     cineq_handle = @(x)f.define_BED_mean_constraint(x,N_obj,nm,rho,BED_mean,n_c(nm),RBE);
                %     options = optimoptions('fmincon','Display','none');
                %     [x,fval,exitflag,output] = fmincon(obj_handle,x,A,b,Aeq,beq,lb,ub,cineq_handle,options);
                %     z_mean{nm} = x;
                % end
            end
        end
        if any(type_constr{nm} == 3)
            for nd = 1:num_DVH(nm)
                if ~isempty(id{nm}{nd})
                    tmp_prod = zeros(numel(id{nm}{nd}),N_obj_c);
                    %tt = zeros(n_c(nm,1));
                    n = 0;
                    for nf = 1:N_obj_c
                        tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                        tmp_prod(:,nf) = tmp_Dij{mod(nf-1,unique_fields)+1,nm}(id{nm}{nd},:)*x0(n+ (1:tmp_index))-lambd_dvh{nm}(id{nm}{nd},nf);
                        %tt = tt+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
                        %z_max{nm}(:,nf) = tmp_prod;
                        n = n+tmp_index;
                    end
                    p = 1 - sqrt(((BED_DVH_max{nm}(nd)/(rho(nm)*N_obj))+(N_obj_c*RBE^2/(2*rho(nm))^2))./sum((tmp_prod+(RBE/(2*rho(nm)))).^2,2));
                    % disp(min(p));
                    % disp(max(p));
                    p = max(0,p);
                    %disp(size(p));
                    z_dvh{nm}(id{nm}{nd},:) = (1-p).*tmp_prod-p.*RBE/(2*rho(nm));
                    tmp = sum(z_dvh{nm},2)*RBE+rho(nm)*sum(z_dvh{nm}.^2,2);
                    % disp(max(tmp));
                    % disp(BED_DVH_max{nm}(nd));
                    % for ii = 1:numel(id{nm}{nd})
                    %     i = id{nm}{nd}(ii);
                    %     tt = 0;
                    %     n = 0;
                    %     for nf = 1:N_obj
                    %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                    %         tt = tt+RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_dvh{nm}(i,nf))+rho(nm)*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_dvh{nm}(i,nf))^2;
                    %         n = n+tmp_index;
                    %     end
                    %     if tt <= BED_DVH_max{nm}(nd)
                    %         n = 0;
                    %         for nf = 1:N_obj
                    %             tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                    %             z_dvh{nm}(i,nf) = tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))-lambd_dvh{nm}(i,nf);
                    %             n = n+tmp_index;
                    %         end
                    %     else
                    %         % Solve the optimization problem
                    %         A = [];
                    %         b = [];
                    %         Aeq = [];
                    %         beq = [];
                    %         ceq = [];
                    %         ub = [];
                    %         lb = [];
                    %         x = ones(N_obj,1);
                    %         obj_handle = @(x)f.define_BED_max_objective(x,N_obj,tmp_Dij,i,nm,lambd,x0);
                    %         cineq_handle = @(x)f.define_BED_max_constraint(x,N_obj,nm,rho,BED_DVH_max{nm}(nd),RBE);
                    %         options = optimoptions('fmincon','Display','off');
                    %         [x,fval,exitflag,output] = fmincon(obj_handle,x,A,b,Aeq,beq,lb,ub,cineq_handle,options);
                    %         z_dvh{nm}(i,:) = x;
                    %     end
                    % end
                end
            end
        end
    end

    %% Update u/x0
    n = 0;
    disp('update x')
    for nf = 1:N_obj_c
        disp(nf)
        tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
        A = ((2*wt_obj(1)/n_target)+(mu_n0*wt_obj(2)/n_target))*Cost_matrix{mod(nf-1,unique_fields)+1}'*Cost_matrix{mod(nf-1,unique_fields)+1} + mu_g*(1/tmp_index)*eye(size(Cost_matrix{mod(nf-1,unique_fields)+1},2));
        b = Cost_matrix{mod(nf-1,unique_fields)+1}'*(2*(wt_obj(1)/n_target)*px*ones(size(Cost_matrix{mod(nf-1,unique_fields)+1},1),1)+mu_n0*(wt_obj(2)/n_target)*(z_n0(:,nf)+lambd_n0(:,nf)));
        if ~isempty(id_min{nf})
            A = A + mu_min*(wt_obj(3)/n_target)*Cost_matrix{mod(nf-1,unique_fields)+1}(id_min{nf},:)'*Cost_matrix{mod(nf-1,unique_fields)+1}(id_min{nf},:);
            b = b + mu_min*(wt_obj(3)/n_target)*Cost_matrix{mod(nf-1,unique_fields)+1}(id_min{nf},:)'*(z_min(id_min{nf},nf)+lambd_min(id_min{nf},nf));
        end
        for m = 1:n_rho
            for ll = 1:length(type_constr{m})
                if type_constr{m}(ll) == 1
                    A = A + mu*(wt_constr{m}(ll)/n_c(m))*tmp_Dij{mod(nf-1,unique_fields)+1,m}'*tmp_Dij{mod(nf-1,unique_fields)+1,m};
                    b = b + mu*(wt_constr{m}(ll)/n_c(m))*tmp_Dij{mod(nf-1,unique_fields)+1,m}'*(z_max{m}(:,nf)+lambd_max{m}(:,nf));
                elseif type_constr{m}(ll) == 2
                    A = A + mu*(wt_constr{m}(ll)/n_c(m))*tmp_Dij{mod(nf-1,unique_fields)+1,m}'*tmp_Dij{mod(nf-1,unique_fields)+1,m};
                    b = b + mu*(wt_constr{m}(ll)/n_c(m))*tmp_Dij{mod(nf-1,unique_fields)+1,m}'*(z_mean{m}(:,nf)+lambd_mean{m}(:,nf));
                end
            end
            if any(type_constr{m} == 3)
                ll = find(type_constr{m} == 3);
                unique_id = [];
                for nd = 1:num_DVH(m)
                    if ~isempty(id{m}{nd})
                        unique_id = [unique_id;id{m}{nd}];
                    end
                end
                unique_id = unique(unique_id);
                if ~isempty(unique_id)
                    A = A + mu_max*(wt_constr{m}(ll)/n_c(m))*tmp_Dij{mod(nf-1,unique_fields)+1,m}(unique_id,:)'*tmp_Dij{mod(nf-1,unique_fields)+1,m}(unique_id,:);
                    b = b + mu_max*(wt_constr{m}(ll)/n_c(m))*tmp_Dij{mod(nf-1,unique_fields)+1,m}(unique_id,:)'*(z_dvh{m}(unique_id,nf)+lambd_dvh{m}(unique_id,nf));
                end
            end
        end
        b = b + mu_g*(1/tmp_index)*(y(n+ (1:tmp_index))+gamma(n+ (1:tmp_index)));
        x0(n+ (1:tmp_index)) = mldivide(A,b);
        n = n + tmp_index;
    end

    %% Update dual variables
    n = 0;
    disp('update dual')
    for nf = 1:N_obj_c
        tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
        lambd_n0(:,nf) = lambd_n0(:,nf) + z_n0(:,nf) - Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index));
        if ~isempty(id_min{nf})
            lambd_min(id_min{nf},nf) = lambd_min(id_min{nf},nf) + z_min(id_min{nf},nf) - Cost_matrix{mod(nf-1,unique_fields)+1}(id_min{nf},:)*x0(n+ (1:tmp_index));
        end
        n = n+tmp_index;
    end
    for m = 1:n_rho
        for ll = 1:length(type_constr{m})
            if type_constr{m}(ll) == 1
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,m},2);
                    lambd_max{m}(:,nf) = lambd_max{m}(:,nf) + z_max{m}(:,nf) - tmp_Dij{mod(nf-1,unique_fields)+1,m}*x0(n+ (1:tmp_index));
                    n = n+tmp_index;
                end
            elseif type_constr{m}(ll) == 2
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,m},2);
                    lambd_mean{m}(:,nf) = lambd_mean{m}(:,nf) + z_mean{m}(:,nf) - tmp_Dij{mod(nf-1,unique_fields)+1,m}*x0(n+ (1:tmp_index));
                    n = n+tmp_index;
                end
            end
        end
        if any(type_constr{m} == 3)
            unique_id = [];
            for nd = 1:num_DVH(m)
                if ~isempty(id{m}{nd})
                    unique_id = [unique_id;id{m}{nd}];
                end
            end
            unique_id = unique(unique_id);
            if ~isempty(unique_id)
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,m},2);
                    lambd_dvh{m}(unique_id,nf) = lambd_dvh{m}(unique_id,nf) + z_dvh{m}(unique_id,nf) - tmp_Dij{mod(nf-1,unique_fields)+1,m}(unique_id,:)*x0(n+ (1:tmp_index));
                    n = n+tmp_index;
                end
            end
        end
    end
    % Update gamma
    gamma = gamma + y - x0;

    %% Update active index set for DVH min constraints
    Viol_DVHmin = 1;
    disp('update ac')
    n = 0;
    for i = 1:N_obj_c
        tmp_index = size(Cost_matrix{mod(i-1,unique_fields)+1},2);
        [id_min{i},Viol] = f.find_active_indices_DVH_min(Cost_matrix{mod(i-1,unique_fields)+1}*x0(n+ (1:tmp_index)),DVH_min_perc,Cost_matrix{mod(i-1,unique_fields)+1},DVH_min);
        Viol_DVHmin = min(Viol_DVHmin,Viol);
        n = n+tmp_index;
    end

    %% Update active index set for DVH max constraints
    Viol_BEDDVH = 0;
    for i = 1:n_rho
        if any(type_constr{i} == 3)
            id{i} = cell(num_DVH(i),1);
            DVH_eval = zeros(n_c(i),1);
            n = 0;
            for nf = 1:N_obj_c
                tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,i},2);
                DVH_eval = DVH_eval + RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index)))+rho(i)*(tmp_Dij{mod(nf-1,unique_fields)+1,i}*x0(n+ (1:tmp_index))).^2;
                n = n+tmp_index;
            end
            DVH_eval = DVH_eval*N_obj;
            if strcmp(ptid,'9306087')
                [id{i},Viol] = f.find_active_indices_DVH_max(DVH_eval,num_DVH(i),DVH_perc{i}/n_c(i),n_c(i),BED_DVH_max{i},id{i});
            else
                [id{i},Viol] = f.find_active_indices_DVH_max(DVH_eval,num_DVH(i),DVH_perc{i},n_c(i),BED_DVH_max{i},id{i});
            end 
            Viol_BEDDVH = max(Viol,Viol_BEDDVH);
        end
    end

    %% Calculate objective function value
    Objective_value = 0;
    disp('calc obj, etc')
    n = 0;
    for i = 1:N_obj_c
        tmp_index = size(Cost_matrix{mod(i-1,unique_fields)+1},2);
        Objective_value = Objective_value+norm(Cost_matrix{mod(i-1,unique_fields)+1}*x0(n+ (1:tmp_index))-px*ones(size(Cost_matrix{mod(i-1,unique_fields)+1},1),1),2)^2;
        n = n+tmp_index;
    end
    Objective_value = Objective_value*N_obj;
    %% Calculate violation of constraints - BED max, BED mean, linear
    Viol_linear = 0;
    n = 0;
    for nf = 1:N_obj_c
        tmp_index = size(Cost_matrix{mod(nf-1,unique_fields)+1},2);
        Viol = max(Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index)),[],"all");
        Viol_linear = max(Viol/px,Viol_linear);
        n = n+tmp_index;
    end
    Viol_linear = max(0,Viol_linear);
    Viol_BEDmax = 0;
    Viol_BEDmean = 0;
    for nm = 1:n_rho
        for ll = 1:length(type_constr{nm})
            if type_constr{nm}(ll) == 1
                tt = zeros(n_c(nm),1);
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                    tmp_prod = tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index));
                    tt = tt+RBE*(tmp_prod)+rho(nm)*(tmp_prod).^2;
                    n = n+tmp_index;
                end
                tt = tt*N_obj;
                Viol = max(tt-BED_max(nm),[],"all");
                Viol_BEDmax = max(Viol_BEDmax,Viol);
                % for i = 1:n_c(nm)
                %     tt = 0;
                %     n = 0;
                %     for nf = 1:N_obj
                %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %         tt = tt+RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index)))+rho(nm)*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index)))^2;
                %         n = n+tmp_index;
                %     end
                %     if tt > BED_max(nm)
                %         Viol_BEDmax = max(Viol_BEDmax,tt-BED_max(nm));
                %     end
                % end
            elseif type_constr{nm}(ll) == 2
                tt = 0;
                n = 0;
                for nf = 1:N_obj_c
                    tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                    tt = tt+sum(RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index)))+rho(nm)*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))).^2);
                    n = n+tmp_index;
                end
                tt = tt*N_obj;
                % for i = 1:n_c(nm)
                %     n = 0;
                %     for nf = 1:N_obj
                %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
                %         tt = tt+RBE*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index)))+rho(nm)*(tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index)))^2;
                %         n = n+tmp_index;
                %     end
                % end
                Viol_BEDmean = max(Viol_BEDmean,tt-n_c(nm)*BED_mean(nm));
                % if tt > n_c(nm)*BED_mean(nm)
                %     Viol_BEDmean = max(Viol_BEDmean,tt-n_c(nm)*BED_mean(nm));
                % end
            end
        end
    end
    disp([k,Objective_value,Viol_linear,Viol_BEDmax,Viol_BEDmean,Viol_BEDDVH,Viol_DVHmin*100]);
    outp.violations{k} = [Objective_value,Viol_linear,Viol_DVHmin,Viol_BEDmax,Viol_BEDmean,Viol_BEDDVH];
end

outp.final_output = [Objective_value,Viol_linear,Viol_DVHmin,Viol_BEDmax,Viol_BEDmean,Viol_BEDDVH];
outp.fluencevector = x0;
fname = strcat('output/BED-ADMM-',string(datetime('now','Format',"yyyy-MM-dd-HH-mm")),'.mat');
if ~exist("output", 'dir')
   mkdir("output")
end
save(fname,"-struct","outp");