% ==============================================================
%%The lines of code here are unique for each test case
% ==============================================================
%% Define data path and parameters
addpath('../utils')
ptid = 'HN02';
folder = ['../' ptid '/'];
f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
t_px = 70; %total prescription dose
nfrac = 35; % number of fraction
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
n_rho = numel(oar)+2;
rho = zeros([n_rho,1]);
rho(1) = cst{15,5}.betaX/cst{15,5}.alphaX;
rho(2) = cst{9,5}.betaX/cst{9,5}.alphaX;
%rho(1) = cst{18,5}.betaX/cst{18,5}.alphaX;
rho(3) = cst{2,5}.betaX/cst{2,5}.alphaX;
rho(4) = cst{11,5}.betaX/cst{11,5}.alphaX;
rho(5) = cst{17,5}.betaX/cst{17,5}.alphaX;
rho(6) = cst{16,5}.betaX/cst{16,5}.alphaX;
rho_target = cst{15,5}.betaX/cst{15,5}.alphaX;
RBE = 1.1;
BED_target = (px*nfrac)*(RBE+rho_target*((px*nfrac)/nfrac));

%% Define arc angles for each fraction
arc_id = cell(nfrac,1); %each cell contains id of the arc associated with that objective
%unique_fields = 6; %There are 6 unique field combinations with 4 fields per fraction as given by arc_id
unique_fields = 1;
for i = 1:nfrac %This can be changed as needed
    %arc_id{i} = [2*(i-1);2*(i+5)]; %;2*i-1;2*(i+5);2*(i+5)+1];
    %arc_id{i} = [6*(i-1)];
    %arc_id{i} = [2*(i-1);2*i-1;2*(i+5);2*(i+5)+1];
    %arc_id{i} = [mod(i-1,unique_fields);mod(i-1,unique_fields)+6;mod(i-1,unique_fields)+12;mod(i-1,unique_fields)+18];
    arc_id{i} = (0:23)';
end
% for i = 1:nfrac
%     disp(arc_id{i}');
% end
K = 30; % number of iterations of the ADMM method

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

% N_obj = 7;
% c_obj = [1; 1; 1; 3; 4; 5; 6; ];
% type_obj = [0;5;8;3;2;1;1;];
% s_obj = [px;1.1*px;px;[30;40;20;20]];
% n_obj = [nan;nan;0.95;0.5;nan;nan;nan;];
% % w_obj = [1;1;100;1;1;1;1];
% w_obj = [1;100;100;1;1;0.1;0.1];
% type_obj = [0;5;8;7;6;5;5;];
% s_obj = [px;1.1*px;px;[30/nfrac;40/nfrac;20/nfrac;20/nfrac]];
% n_obj = [nan;nan;0.95;0.5;nan;nan;nan;];
% w_obj = [1;1;40;1;1;1;1];

N_obj = 7;
c_obj = [1; 1; 1; 3; 4; 5; 6; ];
type_obj = [0;5;8;3;2;2;2;];
s_obj = [px;1.1*px;px;[30;40;50;45;]];
n_obj = [nan;nan;0.95;0.5;nan;nan;nan;];
%w_obj = [1;1;100;1;1;1;1];
w_obj = [1;100;100;1;1;0.1;0.1];

%% Set parameters for the Augmented Lagrangian
mu_n0 = 1e-5;%0.1; %coefficient of penalty for linear constraint: dose to tumor <= 1.1px
mu_min = 1e-5;%0.1; %coefficient of penalty for DVH min constraint for tumor
mu = 1e-3;%0.1; %coefficient of penalty for BED max, BED mean constraint for OAR
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
fname = strcat("Cost_matrix_",int2str(unique_fields),".mat");
if exist(fname,"file")
    disp("Loading data from file...")
    load(fname);
else
    Cost_matrix = cell([min(nfrac,unique_fields),1]);
    num_beams = zeros(min(nfrac,unique_fields),1);
    for i = 1:min(nfrac,unique_fields)
        tmp_indices = [];
        num_beamlets = zeros([numel(arc_id{i}),1]);
        for j = 1:numel(arc_id{i})
            disp('count');
            disp(arc_id{i}(j));
            load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
            num_beamlets(j) = size(dij.physicalDose{1}, 2);
        end
        Cost_matrix{i} = sparse(size(dij.physicalDose{1}, 1),sum(num_beamlets));
        tmp_n = 0;
        for j = 1:numel(arc_id{i})
            load([folder ptid '_' num2str(arc_id{i}(j)*15) '.mat'], 'dij');
            disp('load');
            disp(arc_id{i}(j));
            Cost_matrix{i}(:,tmp_n + (1:num_beamlets(j))) = dij.physicalDose{1}; %Get rows corresponding to tumor voxels, and columns corresponding to arcs
            tmp_n = tmp_n + num_beamlets(j);
        end
        num_beams(i) = sum(num_beamlets);
    end
    save(fname,"Cost_matrix","num_beams");
end

%% Define and initialize the decision variable 'u'. We use x0 in this code.
total_beamlets = 0;
for i = 1:min(nfrac,unique_fields)
    total_beamlets = total_beamlets+size(Cost_matrix{mod(i-1,unique_fields)+1},2);
end
%x0 = ones([total_beamlets,1])*10;
disp('initialize variables');
x0 = zeros([total_beamlets,1]);

%% Initialize all variables in the ADMM method
y = x0;
gamma = y;
z = cell(N_obj,1);
for i = 1:N_obj
    if type_obj(i) == 1 || type_obj(i) == 2 || type_obj(i) == 3 || type_obj(i) == 4
        z{i} = sparse(size(Cost_matrix{1},1),min(nfrac,unique_fields));
        nn = 0;
        for n = 1:min(nfrac,unique_fields)
            tmp_field = mod(n-1,unique_fields)+1;
            tmp_index = num_beams(tmp_field);
            z{i}(c{c_obj(i)},n) = Cost_matrix{tmp_field}(c{c_obj(i)},:)*x0(nn+ (1:tmp_index));
            nn = nn+tmp_index;
        end
    end
end
lambda = z;


%% Define parameters to be saved after the experiment
clear outp
outp.N_obj = N_obj;
outp.c_obj = c_obj;
outp.type_obj = type_obj;
outp.s_obj = s_obj;
outp.n_obj = n_obj;
outp.w_obj = w_obj;
outp.K = K;
outp.mu_names = ["mu_n0","mu_min","mu","mu_max","mu_g"];
outp.mu = [mu_n0,mu_min,mu,mu_max,mu_g];
outp.nfrac = nfrac;
outp.t_px = t_px;
outp.px = px;
outp.arc_id = arc_id;
%outp.viol_names = ["Objective_value","Dmax","DVHmin","Viol_BEDmax","Viol_BEDmean","Viol_BEDDVH"];
outp.violations = cell(K,1);

%% Start the iterations of the ADMM method
for k = 1:K

    %% Find active indices
    disp('update ac')
    id_obj = cell(N_obj,min(nfrac,unique_fields));
    b_obj = cell(N_obj,min(nfrac,unique_fields));
    for i = 1:N_obj
        if type_obj(i) == 0
            for j = 1:min(nfrac,unique_fields)
                id_obj{i,j} = c{c_obj(i)};
                b_obj{i,j} = px*ones(n_c(c_obj(i)),1);
            end
        elseif type_obj(i) == 1 || type_obj(i) == 2 || type_obj(i) == 3 || type_obj(i) == 4
            id = c{c_obj(i)};
            DVH_eval = zeros(n_c(c_obj(i)),min(nfrac,unique_fields));
            n = 0;
            for nf = 1:min(nfrac,unique_fields)
                if mod(nf-1,unique_fields) == 0
                    n = 0;
                end
                tmp_index = num_beams(mod(nf-1,unique_fields)+1);
                tmp_prod = (Cost_matrix{mod(nf-1,unique_fields)+1}(c{c_obj(i)},:)*x0(n+ (1:tmp_index)));
                DVH_eval(:,nf) = DVH_eval(:,nf) + RBE*tmp_prod +rho(c_obj(i))*(tmp_prod).^2;
                n = n+tmp_index;
            end
            DVH_eval1 = zeros(n_c(c_obj(i)),1);
            for nf = 1:nfrac
                tmp_index = mod(nf-1,unique_fields)+1;
                DVH_eval1 = DVH_eval1+DVH_eval(:,tmp_index);
            end
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
                        for j = 1:min(nfrac,unique_fields)
                            id_obj{i,j} =  id(ai:end);
                        end
                    else
                        for j = 1:min(nfrac,unique_fields)
                            id_obj{i,j} = id(ai:id1);
                        end
                    end
                end
                for j = 1:min(nfrac,unique_fields)
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
                        for j = 1:min(nfrac,unique_fields)
                            id_obj{i,j} =  id(1:ai);
                        end
                    else
                        for j = 1:min(nfrac,unique_fields)
                            id_obj{i,j} = id(id1:ai);
                        end
                    end
                end
                for j = 1:min(nfrac,unique_fields)
                    b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                end
            elseif type_obj(i) == 1
                id1 = find(DVH_eval>s_obj(i));
                if ~isempty(id1)
                    for j = 1:min(nfrac,unique_fields)
                        id_obj{i,j} = id(id1);
                        b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                    end
                end
            elseif type_obj(i) == 2
                %disp([i,type_obj(i)]);
                if sum(DVH_eval) > s_obj(i)*n_c(c_obj(i))
                    %disp('here');
                    for j = 1:min(nfrac,unique_fields)
                        id_obj{i,j} = c{c_obj(i)};
                        b_obj{i,j} = z{i}(id_obj{i,j},j)+lambda{i}(id_obj{i,j},j);
                    end
                end
            end
        elseif type_obj(i) == 5 || type_obj(i) == 6 || type_obj(i) == 7 || type_obj(i) == 8
            id = c{c_obj(i)};
            DVH_eval = zeros(n_c(c_obj(i)),min(nfrac,unique_fields));
            n = 0;
            for j = 1:min(nfrac,unique_fields)
                tmp_index = num_beams(mod(j-1,unique_fields)+1);
                DVH_eval = (Cost_matrix{mod(j-1,unique_fields)+1}(c{c_obj(i)},:)*x0(n+ (1:tmp_index)));
                n = n+tmp_index;
                if type_obj(i) == 5
                    id1 = find(DVH_eval>s_obj(i));
                    if ~isempty(id1)
                        id_obj{i,j} = id(id1);
                        b_obj{i,j} = ones(numel(id_obj{i,j}),1)*s_obj(i);
                    end
                elseif type_obj(i) == 6
                    %disp([i,type_obj(i)]);
                    if sum(DVH_eval) > s_obj(i)*n_c(c_obj(i))
                        %disp('here');
                        id_obj{i,j} = c{c_obj(i)};
                        b_obj{i,j} = ones(numel(id_obj{i,j}),1)*s_obj(i);
                    end
                elseif type_obj(i) == 7
                    %disp([i,type_obj(i)]);
                    [eval2, I] = sort(DVH_eval,'descend');
                    id = id(I);
                    ai = ceil(n_obj(i)*n_c(c_obj(i)));
                    if eval2(ai) > s_obj(i)
                        %disp('here');
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
        end
    end



    %% Update y
    disp('update y');
    %y = x0 - gamma;
    y = max(0,x0 - gamma);

    %% Update x0/u
    disp('update x');
    n = 0;
    for j = 1:min(nfrac,unique_fields)
        tmp_index = num_beams(j);
        ATA = zeros(tmp_index,tmp_index);
        b_lin = zeros(tmp_index,1);
        for i = 1:N_obj
            if ~isempty(id_obj{i,j})
                ATA = ATA+(w_obj(i)/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*Cost_matrix{j}(id_obj{i,j},:);
                b_lin=b_lin+(w_obj(i)/n_c(c_obj(i)))*Cost_matrix{j}(id_obj{i,j},:)'*b_obj{i,j};
            end
        end
        ATA = ATA+mu_g*eye(tmp_index);
        b_lin = b_lin+mu_g*(y(n+(1:tmp_index))+gamma(n+(1:tmp_index)));
        x0(n+ (1:tmp_index)) = mldivide(ATA,b_lin);
        n = n + tmp_index;
    end

    %% Update z for active indices
    disp('update z');
    for i=1:N_obj
        if ~isempty(id_obj{i,1})
            if type_obj(i) == 1 || type_obj(i) == 3
                tmp_prod = zeros(numel(id_obj{i,1}),nfrac);
                %tmp_prod = zeros(numel(c{c_obj(i)}),nfrac);
                n = 0;
                for nf = 1:nfrac
                    if nf <= unique_fields
                        tmp_index = num_beams(mod(nf-1,unique_fields)+1);
                        tmp_prod(:,nf) = (Cost_matrix{mod(nf-1,unique_fields)+1}(id_obj{i,1},:)*x0(n+ (1:tmp_index)))- lambda{i}(id_obj{i,1},mod(nf-1,unique_fields)+1);
                        %tmp_prod(:,nf) = (Cost_matrix{mod(nf-1,unique_fields)+1}(c{c_obj(i)},:)*x0(n+ (1:tmp_index)))- lambda{i}(c{c_obj(i)},mod(nf-1,unique_fields)+1);
                        n = n+tmp_index;
                    else
                        tmp_prod(:,nf) = tmp_prod(:,mod(nf-1,unique_fields)+1);
                    end
                end
                p = 1 - sqrt(((s_obj(i)/rho(c_obj(i)))+(nfrac*RBE^2/(2*rho(c_obj(i)))^2))./sum((tmp_prod+(RBE/(2*rho(c_obj(i))))).^2,2));
                p = max(0,p);
                disp(max(p));
                z{i}(id_obj{i,1},:) = (1-p).*tmp_prod(:,1:min(nfrac,unique_fields))-p.*RBE/(2*rho(c_obj(i)));
                %z{i}(c{c_obj(i)},:) = (1-p).*tmp_prod(:,1:min(nfrac,unique_fields))-p.*RBE/(2*rho(c_obj(i)));
            elseif type_obj(i) == 2
                tmp_prod = zeros(numel(id_obj{i,1}),nfrac);
                %tmp_prod = zeros(numel(c{c_obj(i)}),nfrac);
                n = 0;
                for nf = 1:nfrac
                    if nf <= unique_fields
                        tmp_index = num_beams(mod(nf-1,unique_fields)+1);
                        tmp_prod(:,nf) = (Cost_matrix{mod(nf-1,unique_fields)+1}(id_obj{i,1},:)*x0(n+ (1:tmp_index)))- lambda{i}(id_obj{i,1},mod(nf-1,unique_fields)+1);
                        %tmp_prod(:,nf) = (Cost_matrix{mod(nf-1,unique_fields)+1}(c{c_obj(i)},:)*x0(n+ (1:tmp_index)))- lambda{i}(c{c_obj(i)},mod(nf-1,unique_fields)+1);
                        n = n+tmp_index;
                    else
                        tmp_prod(:,nf) = tmp_prod(:,mod(nf-1,unique_fields)+1);
                    end
                end
                p = 1 - sqrt(((s_obj(i)*n_c(c_obj(i))/rho(c_obj(i)))+(n_c(c_obj(i))*nfrac*RBE^2/(2*rho(c_obj(i)))^2))/sum((tmp_prod+(RBE/(2*rho(c_obj(i))))).^2,"all"));
                p = max(0,p);
                z{i}(id_obj{i,1},:) = ((1-p)*tmp_prod(:,1:min(nfrac,unique_fields)))-(p*RBE/(2*rho(c_obj(i))));
                %z{i}(c{c_obj(i)},:) = ((1-p)*tmp_prod(:,1:min(nfrac,unique_fields)))-(p*RBE/(2*rho(c_obj(i))));
            elseif type_obj(i) == 4
                tmp_prod = zeros(numel(id_obj{i,1}),nfrac);
                %tmp_prod = zeros(numel(c{c_obj(i)}),nfrac);
                n = 0;
                for nf = 1:nfrac
                    if nf <= unique_fields
                        tmp_index = num_beams(mod(nf-1,unique_fields)+1);
                        tmp_prod(:,nf) = (Cost_matrix{mod(nf-1,unique_fields)+1}(id_obj{i,1},:)*x0(n+ (1:tmp_index)))- lambda{i}(id_obj{i,1},mod(nf-1,unique_fields)+1);
                        %tmp_prod(:,nf) = (Cost_matrix{mod(nf-1,unique_fields)+1}(c{c_obj(i)},:)*x0(n+ (1:tmp_index)))- lambda{i}(c{c_obj(i)},mod(nf-1,unique_fields)+1);
                        n = n+tmp_index;
                    else
                        tmp_prod(:,nf) = tmp_prod(:,mod(nf-1,unique_fields)+1);
                    end
                end
                p = sqrt(((s_obj(i)/rho(c_obj(i)))+(nfrac*RBE^2/(2*rho(c_obj(i)))^2))./sum((tmp_prod+(RBE/(2*rho(c_obj(i))))).^2,2))-1;
                p = max(0,p);
                z{i}(id_obj{i,1},:) = (1+p).*tmp_prod(:,1:min(nfrac,unique_fields))+p.*RBE/(2*rho(c_obj(i)));
                %z{i}(c{c_obj(i)},:) = (1+p).*tmp_prod(:,1:min(nfrac,unique_fields))+p.*RBE/(2*rho(c_obj(i)));
            end
        end
    end

    %% Update dual variables
    disp('update dual variables');
    gamma = gamma + y - x0;
    % Update lambda{N_obj,1} in R^{nY,6} = lambda + z - Du
    for i = 1:N_obj
        if ~isempty(id_obj{i,1}) && type_obj(i) >=1 && type_obj(i) <= 4
            n = 0;
            for j = 1:min(nfrac,unique_fields)
                tmp_index = num_beams(j);
                lambda{i}(id_obj{i,j},j) = lambda{i}(id_obj{i,j},j)+z{i}(id_obj{i,j},j)-Cost_matrix{j}(id_obj{i,j},:)*x0(n+(1:tmp_index));
                %lambda{i}(c{c_obj(i)},j) = lambda{i}(c{c_obj(i)},j)+z{i}(c{c_obj(i)},j)-Cost_matrix{j}(c{c_obj(i)},:)*x0(n+(1:tmp_index));
                n = n+tmp_index;
            end
        end
    end

    %% Calculate objective function value, and output metrics
    disp('Calculate obj, violations');
    Obj_val = zeros(min(nfrac,unique_fields),N_obj);
    n = 0;
    for j=1:min(nfrac,unique_fields)
        tmp_index = num_beams(j);
        for i=1:N_obj
            if ~isempty(id_obj{i,j})
                Obj_val(j,i) = (0.5*w_obj(i)/n_c(c_obj(i)))*norm((Cost_matrix{j}(id_obj{i,j},:)*x0(n+(1:tmp_index)))-b_obj{i,j},2)^2;
            end
        end
        n = n+tmp_index;
    end
    Total_obj_val = 0;
    %Obj_Violations = zeros(N_obj,1);
    Obj_Violations = sum(Obj_val);
    for nf = 1:nfrac
        tmp_index = mod(nf-1,unique_fields)+1;
        Total_obj_val = Total_obj_val+sum(Obj_val(tmp_index,:));
    end
    if unique_fields == 1
        disp([k Total_obj_val Obj_val]);
        outp.violations{k} = [k Total_obj_val Obj_val];
    else
        disp([k Total_obj_val Obj_Violations]);
        outp.violations{k} = [k Total_obj_val Obj_Violations];
    end
end




outp.final_output = [Total_obj_val Obj_Violations];
outp.fluencevector = x0;


%% Normalize the dose
x0 = max(0,x0);
NF = zeros(min(nfrac,unique_fields),1);
n = 0;
d = zeros(size(Cost_matrix{1},1),min(nfrac,unique_fields));
for nf = 1:min(nfrac,unique_fields)
    [tmp_index1,tmp_index] = size(Cost_matrix{nf});
    y = Cost_matrix{nf}*x0(n+ (1:tmp_index));
    y = y(c{1});
    y2 = sort(y,'descend');
    NF(nf) = px/y2(ceil(n_c(1)*0.95));
    %d_tmp = Cost_matrix{nf}*x0(n+ (1:tmp_index))*NF(nf);
    %Dmax_tmp = max(Dmax_tmp,max(d_tmp(:))/px);
    d(:,nf) = Cost_matrix{mod(nf-1,unique_fields)+1}*x0(n+ (1:tmp_index))*NF(nf);
    n = n+tmp_index;
end
outp.d = d;
outp.NF = NF;

%outp.wt_constr = wt_constr;
fname = strcat('output/BED-ADMM-',string(datetime('now','Format',"yyyy-MM-dd-HH-mm")),'.mat');
if ~exist("output", 'dir')
   mkdir("output")
end
save(fname,"-struct","outp");

