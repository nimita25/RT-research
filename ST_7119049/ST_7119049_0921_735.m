% ==============================================================
%%The lines of code here are unique for each test case
% ==============================================================
%% Define data path and parameters
addpath('../utils')
addpath('../methods')

ptid = '7119049';
folder = ['../' ptid '/'];
f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
K = 30; % number of iterations of the ADMM method
problem = 1; %1: (P1), 2: (P2)

%% Define target and OAR
ctv1 = cst{9,4}{1};   % ptv60, 2*30
body = cst{1,4}{1};
N_oar = 3;
oar = cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20

% Import rho values to calculate BED from physical dose
n_rho = numel(oar)+2;
rho = (1/6)*ones([n_rho,1]);
rho(1) = 1/3;
rho_target = rho(1);

RBE = 1.1;

% Calculate BED for target (rhs of equality constraint)
t_px = 60; %total prescription dose
nfrac = 30; % number of fraction
px = t_px/nfrac; % prescription dose
% Tl = 7;
% Td = 35;
%BED_target = (px*nfrac)*(RBE+rho_target*((px*nfrac)/nfrac)) - max(0,nfrac-Tl-1)*log(2)/Td;
BED_target = (px*nfrac)*(RBE+rho_target*((px*nfrac)/nfrac));

% Set initial T/nfrac; Tl, Td
%nfrac = 45;
nfrac_list = [5 10 15 20 25 30 35 40 45 50];
Tl = 7;
Td = 35;

% Define arc angles for each fraction
arc_id = [0 120 240];

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
%   9 -> maximize target BED objective
%   10 -> minimize OAR BED objective
%   11 -> target BED equality constraint
%   For target: type obj = [0;5;8]
% ==============================================================

% Problem (P1)
if problem == 1
    N_obj = 8;
    c_obj = [1;3;3;3;4;4;5;5;];
    type_obj = [11;6;7;10;7;10;6;10;];
    s_obj = [BED_target;[18;12;0;27;0;20;0;]/nfrac];
    n_obj = [nan;nan;0.3;nan;...
        0.6;nan;nan;nan;];
    %w_obj = [1;1;1;20;1;1;1;20;1;20;1;20;1;20;];
    %w_obj = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;];
    w_obj = [100;1;1;0.05;1;0.05;1;0.05;];
elseif problem == 2
    % Problem (P2)
    N_obj = 5;
    c_obj = [1; 3; 3;  4; 5;];
    type_obj = [9;2;3; 3; 2;];
    s_obj = [0;[18;12;27;20]];
    n_obj = [nan;nan;0.3;0.6;nan;];
    w_obj = [0.001;1;1;1;1;];
end


% N_obj = 13;
% c_obj = [1;1;1;3;3;3;4;4;4;4;5;6;7;];
% type_obj = [0;5;8;3;3;0;3;3;3;0;3;3;3;];
% s_obj = [px;1.1*px;px;[25;35;0;25;35;45;0;25;25;25;]];
% n_obj = [nan;nan;0.95;0.5;0.2;nan;...
%     0.5;0.2;0.1;nan;0.1;0.1;0.5;];
% w_obj = [1;1;1;1;1;0.05;1;1;1;0.05;1;1;1;];



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
%   Cost_matrix - Cost matrix for each unique combo of field; includes beamlets of every
%                 field in the combination
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

% Load Dij matrices
fname = strcat("Cost_matrix.mat");
if exist(fname,"file")
    disp("Loading data from file...")
    load(fname);
else
    Cost_matrix = cell([1,1]);
    tmp_indices = [];
    num_beamlets = zeros([numel(arc_id),1]);
    for j = 1:numel(arc_id)
        disp('count');
        disp(arc_id(j));
        load([folder ptid '_' num2str(arc_id(j)) '.mat'], 'dij');
        num_beamlets(j) = size(dij.physicalDose{1}, 2);
    end
    Cost_matrix{1} = sparse(size(dij.physicalDose{1}, 1),sum(num_beamlets));
    tmp_n = 0;
    for j = 1:numel(arc_id)
        load([folder ptid '_' num2str(arc_id(j)) '.mat'], 'dij');
        disp('load');
        disp(arc_id(j));
        Cost_matrix{1}(:,tmp_n + (1:num_beamlets(j))) = dij.physicalDose{1}; %Get rows corresponding to tumor voxels, and columns corresponding to arcs
        tmp_n = tmp_n + num_beamlets(j);
    end
    num_beams = sum(num_beamlets);
    save(fname,"Cost_matrix","num_beams");
end
% End of parameter setting and loading input

resT = cell(numel(nfrac_list),1);

nfrac_orig = nfrac;

for nfrac = nfrac_list

for i = 1:N_obj
if type_obj(i) == 5 || type_obj(i) == 6 || type_obj(i) == 7 || type_obj(i) == 8
    s_obj(i) = s_obj(i)*(nfrac_orig/nfrac);
end
end

%% Define and initialize the decision variable 'u'. We use x0 in this code.
disp('initialize variables');
x0 = zeros([num_beams,1]);

%% Initialize all variables in the ADMM method
y = x0;
gamma = y;
z = cell(N_obj,1);
for i = 1:N_obj
    if type_obj(i) == 1 || type_obj(i) == 2 || type_obj(i) == 3 || type_obj(i) == 4 || type_obj(i) == 11
        z{i} = sparse(size(Cost_matrix{1},1),1);
        z{i}(c{c_obj(i)},1) = Cost_matrix{1}(c{c_obj(i)},:)*x0;
    end
end
lambda = z;
%nfrac = 25;


%% Define parameters to be saved after the experiment
clear outp
outp.problem = problem;
outp.Tl = Tl;
outp.Td = Td;
outp.N_obj = N_obj;
outp.c_obj = c_obj;
outp.type_obj = type_obj;
outp.s_obj = s_obj;
outp.n_obj = n_obj;
outp.w_obj = w_obj;
outp.K = K;
outp.mu_names = ["mu_n0","mu_min","mu","mu_max","mu_g"];
outp.mu = [mu_n0,mu_min,mu,mu_max,mu_g];
%outp.nfrac = nfrac;
outp.t_px = t_px;
% outp.px = px;
outp.arc_id = arc_id;
outp.BED_target = BED_target;
%outp.viol_names = ["Objective_value","Dmax","DVHmin","Viol_BEDmax","Viol_BEDmean","Viol_BEDDVH"];
outp.violations = cell(K,1);
Best_total_obj_val = 1e10;

%% Start the iterations of the ADMM method
for k = 1:K

    %% Find active indices
    disp('update ac')
    id_obj = cell(N_obj,1);
    b_obj = cell(N_obj,1);
    for i = 1:N_obj
        if type_obj(i) == 0
            for j = 1:1
                id_obj{i,j} = c{c_obj(i)};
                b_obj{i,j} = s_obj(i)*ones(n_c(c_obj(i)),1);
            end
        elseif type_obj(i) == 1 || type_obj(i) == 2 || type_obj(i) == 3 || type_obj(i) == 4
            id = c{c_obj(i)};
            DVH_eval = zeros(n_c(c_obj(i)),1);
            for nf = 1:1
                tmp_prod = (Cost_matrix{1}(c{c_obj(i)},:)*x0);
                DVH_eval(:,nf) = DVH_eval(:,nf) + RBE*tmp_prod +rho(c_obj(i))*(tmp_prod).^2;
            end
            %DVH_eval1 = zeros(n_c(c_obj(i)),1);
            DVH_eval1 = DVH_eval(:,1)*nfrac;
            % for nf = 1:nfrac
            %     DVH_eval1 = DVH_eval1+DVH_eval(:,1);
            % end
            DVH_eval = DVH_eval1;
            if type_obj(i) == 3
                Viol_BEDDVH = 0;
                [eval2, I] = sort(DVH_eval,'descend');
                id = id(I);
                ai = ceil(n_obj(i)*n_c(c_obj(i)));
                if eval2(ai) > s_obj(i)
                    Viol_BEDDVH = max(Viol_BEDDVH, eval2(ai)/s_obj(i));
                    id1 = find(eval2<=0.99*s_obj(i),1,'first');
                    if isempty(id1)
                        for j = 1:1
                            id_obj{i,j} =  id(ai:end);
                        end
                    else
                        for j = 1:1
                            id_obj{i,j} = id(ai:id1);
                        end
                    end
                end
                for j = 1:1
                    b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                end
            elseif type_obj(i) == 4
                Viol_BEDDVH = 0;
                [eval2, I] = sort(DVH_eval,'descend');
                id = id(I);
                ai = ceil(n_obj(i)*n_c(c_obj(i)));
                if eval2(ai) < s_obj(i)
                    Viol_BEDDVH = max(Viol_BEDDVH, eval2(ai)/s_obj(i));
                    id1 = find(eval2<s_obj(i),1,'first');
                    if isempty(id1)
                        for j = 1:1
                            id_obj{i,j} =  id(1:ai);
                        end
                    else
                        for j = 1:1
                            id_obj{i,j} = id(id1:ai);
                        end
                    end
                end
                for j = 1:1
                    b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                end
            elseif type_obj(i) == 1
                id1 = find(DVH_eval>s_obj(i));
                if ~isempty(id1)
                    for j = 1:1
                        id_obj{i,j} = id(id1);
                        b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                    end
                end
            elseif type_obj(i) == 2
                if sum(DVH_eval) > s_obj(i)*n_c(c_obj(i))
                    for j = 1:1
                        id_obj{i,j} = c{c_obj(i)};
                        b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                    end
                end
            end
        elseif type_obj(i) == 5 || type_obj(i) == 6 || type_obj(i) == 7 || type_obj(i) == 8
            id = c{c_obj(i)};
            DVH_eval = zeros(n_c(c_obj(i)),1);
            n = 0;
            for j = 1:1
                tmp_index = num_beams(1);
                DVH_eval = (Cost_matrix{j}(c{c_obj(i)},:)*x0(n+ (1:tmp_index)));
                n = n+tmp_index;
                if type_obj(i) == 5
                    id1 = find(DVH_eval>s_obj(i));
                    if ~isempty(id1)
                        id_obj{i,j} = id(id1);
                        b_obj{i,j} = ones(numel(id_obj{i,j}),1)*s_obj(i);
                    end
                elseif type_obj(i) == 6
                    if sum(DVH_eval) > s_obj(i)*n_c(c_obj(i))
                        id_obj{i,j} = c{c_obj(i)};
                        b_obj{i,j} = ones(numel(id_obj{i,j}),1)*s_obj(i);
                    end
                elseif type_obj(i) == 7
                    [eval2, I] = sort(DVH_eval,'descend');
                    id = id(I);
                    ai = ceil(n_obj(i)*n_c(c_obj(i)));
                    if eval2(ai) > s_obj(i)
                        id1 = find(eval2<=0.99*s_obj(i),1,'first');
                        if isempty(id1)
                            id_obj{i,j} =  id(ai:end);
                        else
                            id_obj{i,j} = id(ai:id1);
                        end
                    end
                    b_obj{i,j} = ones(numel(id_obj{i,j}),1)*s_obj(i);
                elseif type_obj(i) == 8
                    [eval2, I] = sort(DVH_eval,'descend');
                    id = id(I);
                    ai = ceil(n_obj(i)*n_c(c_obj(i)));
                    if eval2(ai) < s_obj(i)
                        id1 = find(eval2<s_obj(i),1,'first');
                        if isempty(id1)
                            id_obj{i,j} =  id(1:ai);
                        else
                            id_obj{i,j} = id(id1:ai);
                        end
                    end
                    b_obj{i,j} = ones(numel(id_obj{i,j}),1)*s_obj(i);
                end
            end
        elseif type_obj(i) == 9 || type_obj(i) == 10
            for j = 1:1
                id_obj{i,j} = c{c_obj(i)};
                b_obj{i,j} = ones(numel(id_obj{i,j}),1);
            end
        elseif type_obj(i) == 11
            id = c{c_obj(i)};
            DVH_eval = zeros(n_c(c_obj(i)),1);
            for nf = 1:1
                tmp_prod = (Cost_matrix{1}(c{c_obj(i)},:)*x0);
                DVH_eval(:,nf) = DVH_eval(:,nf) + RBE*tmp_prod +rho(c_obj(i))*(tmp_prod).^2;
            end
            %DVH_eval1 = zeros(n_c(c_obj(i)),1);
            DVH_eval1 = DVH_eval(:,1)*nfrac - (max(nfrac-1-Tl,0)*log(2)/Td);
            disp([min(DVH_eval1),max(DVH_eval1),mean(DVH_eval1)])
            id1 = find(DVH_eval1< 0.99*s_obj(i) | DVH_eval1> 1.01*s_obj(i));
            if ~isempty(id1)
                for j = 1:1
                    id_obj{i,j} = id(id1);
                    b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                end
            end
            % DVH_eval = sum(DVH_eval1,'all')/n_c(c_obj(i));
            % if DVH_eval < s_obj(i)
            %     id_obj{i,j} = c{c_obj(i)};
            %     b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
            % end
            % DVH_eval = DVH_eval1;
            % id1 = find(DVH_eval ~= s_obj(i));
            % if ~isempty(id1)
            %     for j = 1:1
            %         id_obj{i,j} = id(id1);
            %         b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
            %     end
            % end
        end
    end

    %% Update x0/u
    disp('update x');
    n = 0;
    for j = 1:1
        tmp_index = num_beams(j);
        ATA = zeros(tmp_index,tmp_index);
        b_lin = zeros(tmp_index,1);
        for i = 1:N_obj
            if ~isempty(id_obj{i,j}) && type_obj(i) ~= 9 && type_obj(i) ~= 10
                ATA = ATA+(w_obj(i)/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*Cost_matrix{j}(id_obj{i,j},:);
                b_lin=b_lin+(w_obj(i)/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*b_obj{i,j};
            elseif ~isempty(id_obj{i,j}) && type_obj(i) == 9
                ATA = ATA+((-2*nfrac*rho(c_obj(i))*w_obj(i))/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*Cost_matrix{j}(id_obj{i,j},:);
                b_lin=b_lin+((nfrac*RBE*w_obj(i))/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*b_obj{i,j};
            elseif ~isempty(id_obj{i,j}) && type_obj(i) == 10
                ATA = ATA+((2*nfrac*rho(c_obj(i))*w_obj(i))/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*Cost_matrix{j}(id_obj{i,j},:);
                b_lin=b_lin+((-nfrac*RBE*w_obj(i))/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*b_obj{i,j};
            end
        end
        ATA = ATA+mu_g*eye(tmp_index);
        b_lin = b_lin+mu_g*(y(n+(1:tmp_index))+gamma(n+(1:tmp_index)));
        x0(n+ (1:tmp_index)) = mldivide(ATA,b_lin);
        n = n + tmp_index;
    end

    %% Update y
    disp('update y');
    %y = x0 - gamma;
    y = max(0,x0 - gamma);

    %% Update z for active indices
    disp('update z');
    for i=1:N_obj
        if ~isempty(id_obj{i,1})
            if type_obj(i) == 1 || type_obj(i) == 3
                tmp_prod = zeros(numel(id_obj{i,1}),1);
                tmp_prod(:,1) = (Cost_matrix{1}(id_obj{i,1},:)*x0)- lambda{i}(id_obj{i,1},1);
                % for nf = 2:nfrac
                %     tmp_prod(:,nf) = tmp_prod(:,1);
                % end
                p = 1 - sqrt(((s_obj(i)/rho(c_obj(i)))+(nfrac*RBE^2/(2*rho(c_obj(i)))^2))./(nfrac*(tmp_prod+(RBE/(2*rho(c_obj(i))))).^2));
                p = max(0,p);
                z{i}(id_obj{i,1},:) = (1-p).*tmp_prod(:,1)-p.*RBE/(2*rho(c_obj(i)));
            elseif type_obj(i) == 2
                tmp_prod = zeros(numel(id_obj{i,1}),1);
                tmp_prod(:,1) = (Cost_matrix{1}(id_obj{i,1},:)*x0)- lambda{i}(id_obj{i,1},1);
                % for nf = 2:nfrac
                %     tmp_prod(:,nf) = tmp_prod(:,1);
                % end
                p = 1 - sqrt(((s_obj(i)*n_c(c_obj(i))/rho(c_obj(i)))+(n_c(c_obj(i))*nfrac*RBE^2/(2*rho(c_obj(i)))^2))/sum((nfrac*(tmp_prod+(RBE/(2*rho(c_obj(i))))).^2)),'all');
                p = max(0,p);
                z{i}(id_obj{i,1},:) = ((1-p)*tmp_prod(:,1))-(p*RBE/(2*rho(c_obj(i))));
            elseif type_obj(i) == 4
                tmp_prod = zeros(numel(id_obj{i,1}),1);
                tmp_prod(:,1) = (Cost_matrix{1}(id_obj{i,1},:)*x0)- lambda{i}(id_obj{i,1},1);
                p = sqrt(((s_obj(i)/rho(c_obj(i)))+(nfrac*RBE^2/(2*rho(c_obj(i)))^2))./(nfrac*(tmp_prod+(RBE/(2*rho(c_obj(i))))).^2))-1;
                p = max(0,p);
                z{i}(id_obj{i,1},:) = (1+p).*tmp_prod(:,1)+p.*RBE/(2*rho(c_obj(i)));
            elseif type_obj(i) == 11
                % tmp_prod = zeros(numel(id_obj{i,1}),1);
                % tmp_prod(:,1) = (Cost_matrix{1}(id_obj{i,1},:)*x0)- lambda{i}(id_obj{i,1},1);
                tmp_prod = zeros(numel(c{c_obj(i)}),1);
                tmp_prod(:,1) = (Cost_matrix{1}(c{c_obj(i)},:)*x0)- lambda{i}(c{c_obj(i)},1);
                rhs = (((BED_target+(max(nfrac-1-Tl,0)*log(2)/Td))/(rho(c_obj(i))*nfrac))+(RBE/(2*rho(c_obj(i))))^2)*n_c(c_obj(i));
                centr = RBE/(2*rho(c_obj(i)));
                p = sqrt(rhs/sum((tmp_prod+centr).^2,'all')) -1;
                %p = max(0,p);
                % z{i}(id_obj{i,1},:) = (1+p)*tmp_prod(:,1) + p*centr;
                % disp([max(z{i}(id_obj{i,1},:)),min(z{i}(id_obj{i,1},:))])
                z{i}(c{c_obj(i)},:) = (1+p)*tmp_prod(:,1) + p*centr;
                disp([max(z{i}(c{c_obj(i)},:)),min(z{i}(c{c_obj(i)},:))])
                %z{i}(id_obj{i,1},:) = ones(numel(id_obj{i,1}),1)*(sqrt(rhs) - centr);
            end
        end
    end

    %% Update T/nfrac
    if problem == 1 && k > 50
        %tmp_prod = Cost_matrix{1}(c{c_obj(1)},:)*x0;
        tmp_prod = z{1}(c{c_obj(1)},:);
        sl = (1/n_c(c_obj(1)))*(sum((RBE*tmp_prod)+(rho(c_obj(1))*tmp_prod.^2),'all'));
        disp(nfrac*(sl-(log(2)/Td)))
        disp(BED_target - (Tl+1)*log(2)/Td)
        disp(sl - log(2)/Td)
        if nfrac <= Tl+1
            if nfrac*sl > BED_target
                disp('1')
                nfrac = BED_target/sl;
            end
        else
            if nfrac*(sl-(log(2)/Td)) == BED_target - (Tl+1)*log(2)/Td && sl - log(2)/Td <= 0
                disp('2')
                nfrac = BED_target/sl;
            elseif nfrac*(sl-(log(2)/Td)) > BED_target - (Tl+1)*log(2)/Td && sl - log(2)/Td <= 0
                disp('3')
                nfrac = BED_target/sl;
            elseif nfrac*(sl-(log(2)/Td)) > BED_target - (Tl+1)*log(2)/Td && sl - log(2)/Td > 0
                tmpBED = sl*(Tl+1);
                if BED_target <= tmpBED
                    disp('4')
                    nfrac = BED_target/sl;
                else
                    disp('5')
                    nfrac = (BED_target - (Tl+1)*log(2)/Td)/(sl-(log(2)/Td));
                end
            end
        end
    elseif problem == 2
        if k > 50
            tmp_prod = Cost_matrix{1}(c{c_obj(1)},:)*x0;
            %tmp_prod = z{1}(c{c_obj(1)},:);
            sl = (1/n_c(c_obj(1)))*(sum((RBE*tmp_prod)+(rho(c_obj(1))*tmp_prod.^2),'all'));
            ine = 0;
            for i = 1:N_obj
                if type_obj(i) == 1 || type_obj(i) == 2 || type_obj(i) == 3
                    if ~isempty(id_obj{i,1})
                        ine = 1;
                        break
                    end
                end
            end
            disp(sl - log(2)/Td)
            disp(ine)
            if ine == 1 && sl - log(2)/Td <= 0
                nfrac = min(nfrac,Tl+1);
            elseif ine == 0
                tmpT = 1e10;
                for i = 1:N_obj
                    if type_obj(i) == 1 
                        lhs =z{i}(c{c_obj(i)},:);
                        tmpT = min(tmpT,min(s_obj(i)./max(0,((RBE*lhs)+(rho(c_obj(i))*lhs.^2)))));
                    elseif type_obj(i) == 2
                        lhs = z{i}(c{c_obj(i)},:);% Cost_matrix{1}(c{c_obj(i)},:)*x0;
                        tmpT = min(tmpT,(s_obj(i)*n_c(c_obj(i)))/max(0,(sum((RBE*lhs)+(rho(c_obj(i))*lhs.^2),'all'))));
                    elseif type_obj(i) == 3
                        lhs = z{i}(c{c_obj(i)},:);%Cost_matrix{1}(c{c_obj(i)},:)*x0;
                        lhs = (RBE*lhs)+(rho(c_obj(i))*lhs.^2);
                        id = c{c_obj(i)};
                        [eval2, I] = sort(lhs,'descend');
                        ai = ceil(n_obj(i)*n_c(c_obj(i)));
                        tmpT = min(tmpT,min(s_obj(i)./max(0,((RBE*eval2(ai:end))+(rho(c_obj(i))*eval2(ai:end).^2)))));
                        %tmpid = id(ai:end);
                    end
                end
                if sl - log(2)/Td > 0
                    nfrac = min(45,tmpT);
                else
                    nfrac = min(tmpT,Tl+1);
                end
            end
        end
    end

    %% Update dual variables
    disp('update dual variables');
    gamma = gamma + y - x0;
    % Update lambda{N_obj,1} in R^{nY,6} = lambda + z - Du
    for i = 1:N_obj
        if ~isempty(id_obj{i,1}) && (type_obj(i) ==1 || type_obj(i) == 2 || type_obj(i) == 3 || type_obj(i) == 4 || type_obj(i) == 11)
            n = 0;
            for j = 1:1
                tmp_index = num_beams(j);
                %lambda{i}(id_obj{i,j},j) = lambda{i}(id_obj{i,j},j)+z{i}(id_obj{i,j},j)-Cost_matrix{j}(id_obj{i,j},:)*x0(n+(1:tmp_index));
                lambda{i}(c{c_obj(i)},j) = lambda{i}(c{c_obj(i)},j)+z{i}(c{c_obj(i)},j)-Cost_matrix{j}(c{c_obj(i)},:)*x0(n+(1:tmp_index));
                n = n+tmp_index;
            end
        end
    end

    %% Calculate objective function value, and output metrics
    disp('Calculate obj, violations');
    Obj_val = zeros(1,N_obj);
    n = 0;
    for j=1:1
        tmp_index = num_beams(j);
        for i=1:N_obj
            if ~isempty(id_obj{i,j}) && type_obj(i) ~= 9 && type_obj(i) ~= 10
                Obj_val(j,i) = (0.5*w_obj(i)/n_c(c_obj(i)))*norm((Cost_matrix{j}(id_obj{i,j},:)*x0(n+(1:tmp_index)))-b_obj{i,j},2)^2;
            elseif ~isempty(id_obj{i,j}) && type_obj(i) == 9
                tmp_prod = Cost_matrix{j}(id_obj{i,j},:)*x0(n+(1:tmp_index));
                Obj_val(j,i) = ((-nfrac*w_obj(i))/n_c(c_obj(i)))*sum(RBE*tmp_prod + rho(c_obj(i))*tmp_prod.^2,'all')+ w_obj(i)*(max(nfrac-1-Tl,0)*log(2)/Td);
            elseif ~isempty(id_obj{i,j}) && type_obj(i) == 10
                tmp_prod = Cost_matrix{j}(id_obj{i,j},:)*x0(n+(1:tmp_index));
                Obj_val(j,i) = ((nfrac*w_obj(i))/n_c(c_obj(i)))*sum(RBE*tmp_prod + rho(c_obj(i))*tmp_prod.^2,'all');
            end
        end
        n = n+tmp_index;
    end
    Total_obj_val = sum(Obj_val);
    disp([k Total_obj_val Obj_val nfrac]);
    outp.violations{k} = [k Total_obj_val Obj_val nfrac];
    if Total_obj_val < Best_total_obj_val && k > 25
        if (problem == 1 && Total_obj_val > 0) || problem == 2
            Best_total_obj_val = Total_obj_val;
            outp.best_output = [k Total_obj_val Obj_val];
            outp.best_fluencevector = x0;
            outp.best_nfrac = nfrac;
        end
    end

end


outp.final_output = [Total_obj_val Obj_val];
outp.fluencevector = x0;
outp.nfrac = nfrac;

% Normalize the dose
x0 = max(0,x0);
d = Cost_matrix{1}*x0;
BED2 = zeros(size(d,1),1);
BED2 = nfrac*(RBE*d + rho(1)*d.^2);
i=1;
BED2(c{i}) = BED2(c{i}) - max(nfrac-Tl-1,0)*log(2)/Td;
% for i = 1:numel(c)
%     BED2(c{i}) = RBE*d(c{i}) + rho(i)*d(c{i}).^2;
%     BED2(c{i}) = BED2(c{i})*nfrac;
%     if i == 1
%         BED2(c{i}) = BED2(c{i}) - max(nfrac-Tl-1,0)*log(2)/Td;
%     end
% end

meanBEDtarget = mean(BED2(c{1}));
px = sqrt((RBE/rho(1))^2+(4*(meanBEDtarget+(max(nfrac-Tl-1,0)*log(2)/Td)))/(nfrac*rho(1)))-(RBE/rho(1));
px = px/2;


[D95, Dmax, CI, BEDmean_lung, BED30_lung, BEDmean_heart, BED60_heart, ...
    BEDmean_eso] =  calcparaBED_7119049(d, BED2, px, c);

resd =[D95, Dmax, CI, BEDmean_lung, BED30_lung, BEDmean_heart, BED60_heart, ...
    BEDmean_eso];

outp.output_meanBEDtarget = meanBEDtarget;
outp.d = d;
outp.resd = resd;
outp.output_px = px;

idx = find(nfrac_list==nfrac);
resT{idx} = outp;

end

%outp.wt_constr = wt_constr;
fname = strcat('output/ST-',string(datetime('now','Format',"yyyy-MM-dd-HH-mm")),'.mat');
if ~exist("output", 'dir')
   mkdir("output")
end
%save(fname,"-struct","outp");
save(fname,"resT");

