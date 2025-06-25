
%clc;clear;close all;
ptid = '7119049';
folder = ['../' ptid '/'];
px = 2; % prescription dose
nfrac = 30; % number of fraction
mu_min = 5; % MMU threshold



%% Define plan objective
load([folder ptid '.mat'],'ct','cst');
% Extract CTV and body information
ctv1 = cst{9,4}{1};   % ptv60, 2*30
body = cst{1,4}{1};
N_oar = 3;
oar = cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20

load([folder ptid '.mat'],'ct','cst');
ctv = cell(1,1);
ctv{1} = ctv1;
n_oar = zeros(N_oar,1);

N_Dij = 1;
N_obj = N_Dij;
chioce_case = 2;
%NNZ = 40; % number of nonzeros

% %% Load Data
% Dij = cell(N_Dij,1);
% x_init = [];
% for k = 1:N_Dij
%     load(['./results/' ptid '_ADMM_mmu_' num2str(mu_min) '.mat']);
%     x_init = x0;
%     load([folder ptid '_v' num2str(k) '.mat']);
% 
%     [nY,nX] = size(Dij0);
%     Dij{k} = Dij0;
% end
% [nY,nX] = size(Dij0);


%% Define  delivery angles
if 1
    id = [0 120 240];
    N_iter = 3;
    N_iter_admm = 50;
else
    id = 0:15:345;
    N_iter = 50;
end

%% Load influence matrix
load([folder ptid '_' num2str(id(1)) '.mat'], 'dij', 'stf');
m = size(dij.physicalDose{1}, 1);
N = zeros(numel(id), 1);
N(1) = size(dij.physicalDose{1}, 2);
for i = 2:numel(id)
    load([folder ptid '_' num2str(id(i)) '.mat'], 'dij');
    N(i) = size(dij.physicalDose{1}, 2);
end

Dij = sparse(m, sum(N));
n = 0;
for i = 1:numel(id)
    load([folder ptid '_' num2str(id(i)) '.mat'], 'dij' , 'stf');
    Dij(:, n + (1:N(i))) = dij.physicalDose{1};
    n = n + N(i);
end
[nY, nX] = size(Dij);
Dij = {Dij};

%% 2. Find energy beam intensity correspondence (not necessary, may be used for mixed integer)
%load(['./results/' ptid '_ADMM_mmu_' num2str(mu_min) '.mat']);
nf = numel(id);%number of angle
pe = cell(nf,1);
eg = cell(nf,1);
pe_unique = cell(nf,1);
ii = 0;
for i = 1:nf
    load([folder ptid '_' num2str(id(i)) '.mat'], 'dij' , 'stf');
    pe{i} = zeros([stf.totalNumOfBixels 1]);
    nr = numel(stf.numOfBixelsPerRay);
    ie = 0;
    for j = 1:nr
        ne = stf.numOfBixelsPerRay(j);
        pe{i}(ie+(1:ne)) = stf.ray(j).energy;
        ie = ie + ne;
    end
    pe_unique{i} = unique(pe{i});
    eg(i) = {cell([numel(pe_unique{i}),1])};
    for j = 1:numel(pe_unique{i})
        eg{i}{j} = find(pe{i} == pe_unique{i}(j)) + ii;
    end
    ii = ii + stf.totalNumOfBixels;
end
ne_unique = 0;
pe2 = [];
for i = 1:nf
    ne_unique = ne_unique+numel(pe_unique{i});
    pe2 = [pe2;pe{i}];
end
id_gs = zeros([nX 1]);
n_gs = zeros([ne_unique 1]);
ii = 0; jj = 0;
for i = 1:nf
    for j = 1:numel(pe_unique{i})
        ii = ii + 1;
        n_gs(ii) = numel(eg{i}{j});
        id_gs(jj + (1:n_gs(ii))) = eg{i}{j};
        jj = jj + n_gs(ii);
    end
end
for i = 1:N_Dij % N_dij : the number of the scenarios
    Dij{i} = Dij{i}(:,id_gs);
end

for NNZ = 40:5:50
    mu_list = [1e-11 1e-9 1e-5 1e-1 1e1];
    for mm = mu_list
    [nY,nX] = size(Dij{1});
    
    %%  Define optimization objectives
    %% 
    % ==============================================================
    %   c     - Row index for different structures.
    %   N_c   - Number of the structures.
    %   n_c   - Array store the number of voxels for each structure.
    % ==============================================================
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
    
    % ==============================================================
    %   N_obj     - Number of objectives.
    %   w_obj     - Objective weight.
    %   s_obj     - Prescription dose.
    %   n_obj     - Array stores the number of active index for each objective.
    %   c_obj     - Strucuture index for each objective.
    %   id_obj    - Active index for each DVH objective.
    %   type_obj  - DVH type of each objective.
    % ==============================================================
    
    
    px0 = px * nfrac;
    N_obj = 12;
    type_obj = [0;2;3;0;3;1;1;0;1;1;0;0];
    w_obj = [1;1;1;0.1;1;1;1;0.005;1;1;0.01;0.01]; 
    s_obj = [px;px;px*1.1;0;px;[18;12;0;27;20;0;0]/nfrac];
    n_obj = round([nan;n_c(1)*0.95;nan;nan;nan;n_c(3)*0.5;n_c(3)*0.3;nan;n_c(4)*0.6;n_c(5)*0.5;nan;nan;]);
    c_obj = [1;1;1;2;2;3;3;3;4;5;4;5];
    id_obj = cell(N_obj, 1);
    % AmX = @AmX_robust;
    % AtmX = @AtmX_robust;
    % Update_ac = @update_ac_robust;
    % Calc_obj_dvh = @calc_obj_dvh_robust;
    AmX = @AmX_v3; 
    AtmX = @AtmX_v3;
    Calc_obj_dvh = @calc_obj_dvh;
    Update_ac =  @update_ac;
    var_plot = struct('n_c', n_c, 'dmax', px * 1.2);
    para = struct('Dij', {Dij}, 'nX', nX, 'nY', nY, 'N_c', N_c, 'n_c', n_c, 'c', {c}, 'isC', uint32(0),...
        'N_obj', N_obj, 'type_obj', type_obj, 'w_obj', w_obj, 's_obj', s_obj, 'n_obj', n_obj, 'id_obj', {id_obj}, 'c_obj', c_obj);
    %ip = struct('N_iter', [], 'nplot', 10000, 'var_plot', var_plot, 'isC', uint32(0));
    
    
    %% pre-optimization 
    N_obj = 2;
    type_obj = [0; 0;];
    w_obj = [1; 0.1;];
    s_obj = [px; 0;];
    n_obj = round([nan; nan;]);
    c_obj = [1; 2;];
    id_obj = cell(N_obj, 1);
    id_obj(1) = c(1);
    id_obj(2) = c(2);
    para0 = struct('Dij', {Dij}, 'nX', nX, 'nY', nY, 'N_c', N_c, 'n_c', n_c, 'c', {c}, 'isC', uint32(0),...
        'N_obj', N_obj, 'type_obj', type_obj, 'w_obj', w_obj, 's_obj', s_obj, 'n_obj', n_obj, 'id_obj', {id_obj}, 'c_obj', c_obj);
    x0 = ones([nX 1], 'single');
    Ni = 100;
    ns = numel(n_gs);
    wp = struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(n_gs),'isC',uint32(0),'w',ones([nX 1],'single'));
    ip = struct('Ni',Ni,'N_iter',N_iter,'ns',ns,'update_ac',Update_ac,'calc_obj_dvh',Calc_obj_dvh,'isC',uint32(0),'nplot',10000,'var_plot', var_plot);
    x0 = ones([nX 1],'single');
    maxAtA1 = norm(AtmX(AmX(x0,para0),para0))/norm(x0);
    
    %% Optimization
    var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
        'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', [],'wpr',wp);
    
    nnz_x = 10;
    ip.mup = maxAtA1*1e-2;%1e-1
    ip.mu = mm;%1e-9;%maxAtA2*1e-1;
    ip.mu_min = mu_min;
    ip.nnz_x = nnz_x;
    % load(['./Results_' ptid '/res_' ptid  '.mat']);
    % load(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_CARD_50.mat']);
    % x_init = x0;
    %load(['.\Results_' ptid '\res2308_' ptid '_NNZ_' int2str(NNZ) '_CARD_50.mat'],'sg');
    x_init = zeros([nX 1], 'single');
    %x_init(1:NNZ) = 1;
    tic; [x0,s,sg] = admm_ip_2308(ip,var_CG,NNZ,x_init); toc;

    %% Call ADMM method to optimize over active energy layers
    nn = 0;
    %var_CG.var_AtA.Dij = [];
    tmp_Dij = [];
    for ii = 1:ns
        if ismember(ii,sg)
            tmp_Dij = [tmp_Dij Dij{1}(:,nn+(1:n_gs(ii)))];
        end
        nn = nn+n_gs(ii);
    end
    tmp_Dij = {tmp_Dij};
    [nY, nX] = size(tmp_Dij{1});
    x_init = zeros([nX 1], 'single');

    para.Dij = tmp_Dij;
    para.nX = nX;
    para.nY = nY;
    para0.Dij = tmp_Dij;
    para0.nX = nX;
    para0.nY = nY;
    
    %% Optimization
    % mup = 0.01;
    % ip.mup = maxAtA * mup;
    ip.N_iter = N_iter_admm;
    var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
        'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', []);
    x1 = admm_mmu_withDVH(ip, var_CG);

    %Scale up x1 to original size
    nX = size(Dij{1},2);
    x0 = zeros(nX,1);
    nn = 0;
    nn1 = 0;
    for ii = 1:ns
        if ismember(ii,sg)
            x0(nn+(1:n_gs(ii))) = x1(nn1+(1:n_gs(ii)));
            nn1 = nn1+n_gs(ii);
        end
        nn = nn+n_gs(ii);
    end
    
    % Project to MMU 
    xp = double(x0);
    xp(xp < mu_min / 2) = 0;
    xp(intersect(find(xp >= mu_min / 2), find(xp < mu_min))) = mu_min;
    
    %% Calculate plan parameters
    % plan normalization %
    n_ctv95 = round(n_c(1)*0.95);
    if 1
        y0 = Dij{1}*double(xp);
        y = y0(c{1});
        y2 = sort(y,'descend');
        factor = px/y2(n_ctv95);
        x = xp*factor;
    else
        x = xp;
    end
    
    d = Dij{1}*double(x);
    y = d(c{1});
    y2 = sort(y,'descend');
    D95 = y2(n_ctv95)/px;
    Dmax = y2(1)/px;
    d3d = reshape(d,ct.cubeDim);

    para.Dij = Dij;
    [nY, nX] = size(Dij{1});
    para.nY = nY;
    para.nX = nX;
    [tmp,v] = Update_ac(x,para);
    [obj,D] = Calc_obj_dvh(x,v);
    [mean(sum(obj)) D95 Dmax]
    %dlmwrite('GridSearch_SmartInit.csv', [i_mup i_mu mean(sum(obj)) D95 Dmax],'delimiter',',','-append');
    obj_total = mean(sum(obj));

    [D95, Dmax, CI, Dmean_lung, V20_lung, Dmean_heart, V30_heart, Dmean_esoph,Dmean_body] = calcpara_7119049_2308(d, px, c);
    resd = [D95, Dmax, CI, Dmean_lung, V20_lung, Dmean_heart, V30_heart, Dmean_esoph,Dmean_body];

    % Save output
    clear outp
    outp.obj_total = obj_total;
    outp.x = x;
    outp.D95 = D95;
    outp.Dmax = Dmax;
    outp.d = d;
    outp.d3d = d3d;
    outp.n_gs = n_gs;
    outp.sg = sg;
    outp.obj = obj;
    outp.x0 = x0;
    outp.xp = xp;
    outp.s = s;
    outp.resd = resd;
    outp.mu = ip.mu;
    outp.mu_min = mu_min;
    %fname =strcat('.\Results_', ptid, '\res_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter_admm), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
    fname =strcat('.\Results_', ptid, '\res2908_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter_admm), string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
    %fname =strcat('.\Results_', ptid, '\res_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter_admm), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
    %fname =strcat('.\Results_', ptid, '\res2908_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter_admm), '.mat');
    %save(fname,'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
    save(fname,"-struct","outp");


    % save(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_MIP_' int2str(N_iter) '.mat'],'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
    %save(['.\Results_' ptid '\res_' ptid '_MIP_50_10.mat'],'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
    % if N_Dij == 1
    %     save(['.\Results_' ptid '\res_' ptid '_' num2str(i_mup) '_' num2str(i_mu) '_MIP.mat'],'x','x0','obj','d','Dmax','n_gs');
    % end
    end
end