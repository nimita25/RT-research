for NNZ = 5:5:35
    %clc;clear;close all;
    ptid = '51477044';
    folder = ['../' ptid '/'];
    px = 1.8; % prescription dose
    nfrac = 25; % number of fraction
    mu_min = 5; % MMU threshold
    
    
    
    %% Define plan objective
    load([folder ptid '.mat'],'ct','cst');
    % Extract CTV and body information
    ctv1 = cst{10,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
    body = cst{18,4}{1};
    N_oar = 5;
    oar = cell(N_oar,1);
    oar(1) = cst{9,4}(1); % bladder D50: 25Gy; D20: 35Gy.
    oar(2) = cst{17,4}(1); % rectum D50: 25Gy; D20: 35Gy; D10: 45Gy.
    oar(3) = cst{11,4}(1); % femhead_lt D10: 25Gy.
    oar(4) = cst{12,4}(1); % femhead_rt D10: 25Gy.
    oar(5) = cst{14,4}(1); % penilebulb D50: 25Gy.
    
    load([folder ptid '.mat'],'ct','cst');
    ctv = cell(1,1);
    ctv{1} = ctv1;
    n_oar = zeros(N_oar,1);
    
    N_Dij = 1;
    N_obj = N_Dij;
    choice_case = 2;
    %NNZ = 30; % number of nonzeros
    
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
        id = [90 270];
        N_iter = 50;
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
    N_obj = 15;
    w_obj = [1;1;1;0.1;1;1;1;0.05;1;1;1;0.05;1;1;1;];
    s_obj = [45;45;45*1.1;0;45;[25;35;0;25;35;45;0;25;25;25;]]*px/45;
    n_obj = round([nan;n_c(1)*0.98;nan;nan;nan;n_c(3)*0.5;n_c(3)*0.2;nan;...
        n_c(4)*0.5;n_c(4)*0.2;n_c(4)*0.1;nan;n_c(5)*0.1;n_c(6)*0.1;n_c(7)*0.5;]);
    c_obj = [1;1;1;2;2;3;3;3;4;4;4;4;5;6;7;];
    id_obj = cell(N_obj, 1);
    type_obj = [0;2;3;0;3;1;1;0;1;1;1;0;1;1;1;];
    AmX = @AmX_robust;
    AtmX = @AtmX_robust;
    Update_ac = @update_ac_robust;
    Calc_obj_dvh = @calc_obj_dvh_robust;
    % Update_ac = @update_ac;
    % Calc_obj_dvh = @calc_obj_dvh;
    var_plot = struct('n_c', n_c, 'dmax', px * 1.2);
    para = struct('Dij',{Dij},'nX',nX,'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),'N_Dij',N_Dij,...
        'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,'s_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'c_obj',c_obj);
    
    
    %% pre-optimization 
    N_obj = 2;
    type_obj = [0; 0;];
    w_obj = [1;0.1;];
    w_obj = repmat(w_obj,1,N_obj);
    s_obj = [px;0;];
    n_obj = round([nan;nan;]);
    c_obj = [1;2;];
    id_obj = cell(N_obj,N_Dij);
    for j = 1:N_Dij
        for i = 1:N_obj
            id_obj{i,j} = c{c_obj(i)};
        end
    end
    para0.id_obj = id_obj;
    para0 = para;
    para0.w_obj = w_obj;
    para0.type_obj = type_obj ;
    para0.s_obj  = s_obj ;
    para0.n_obj = n_obj;
    para0.c_obj = c_obj;
    para0.id_obj = id_obj;
    para0.N_obj = N_obj;
    Ni = 100;
    var_gs = struct('Nx',uint32(nX),'na',uint32(numel(n_gs)),'nx',uint32(n_gs));
    ns = numel(n_gs);
    wp = struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(n_gs),'isC',uint32(0),'w',ones([nX 1],'single'));
    wps = struct('isC',uint32(0));
    var_CG = struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(1),'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
        'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,'wps',wps,'mu_xs',[],'Wxr',@GX,'Wtxr',@GtX,'wpr',wp,'mu_xr',[]);
    ip = struct('Ni',Ni,'N_iter',N_iter,'ns',ns,'update_ac',Update_ac,'calc_obj_dvh',Calc_obj_dvh,'isC',uint32(0),'nplot',10000);
    x0 = ones([nX 1],'single');
    maxAtA1 = norm(AtmX(AmX(x0,para0),para0))/norm(x0);
    maxAtA2 = norm(AtmX(AmX(x0,para0),para0))/norm(GtX(GX(x0,wp),wp));
    
    
    %% Optimization
    nnz_x = 10;
    ip.mup = maxAtA1*1e-2;%1e-1
    ip.mu = maxAtA2*1e-1;
    ip.mu_min = mu_min;
    ip.nnz_x = nnz_x;
    x_init = zeros([nX 1], 'single');
    % load(['.\Results_' ptid '\res_' ptid '_ADMM_60.mat']);
    % x_init = x0;
    tic; [x0,sg] = admm_mmu_gs(ip,var_CG,NNZ,x_init); toc;
    
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

    %Calculate output metrics    
    d = Dij{1}*double(x);
    y = d(c{1});
    y2 = sort(y,'descend');
    D95 = y2(n_ctv95)/px;
    Dmax = y2(1)/px;
    d3d = reshape(d,ct.cubeDim);
    
    
    [tmp,v] = Update_ac(x,para);
    [obj,D] = Calc_obj_dvh(x,v);
    [mean(sum(obj)) D95 Dmax]
    obj_total = mean(sum(obj));

    [D95, Dmax, CI, Dmean_oar1, V10_oar1, Dmean_oar2, V10_oar2, Dmean_oar3, Dmean_oar5, Dmean_body] = calcpara_51477044_2308(d, px, c);
    resd = [D95, Dmax, CI, Dmean_oar1, V10_oar1, Dmean_oar2, V10_oar2, Dmean_oar3, Dmean_oar5, Dmean_body];

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
    outp.resd = resd;
    fname =strcat('.\Results_', ptid, '\res2308_', ptid, '_NNZ_', int2str(NNZ), '_CARD_', int2str(N_iter), '.mat');
    %save(fname,'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
    save(fname,"-struct","outp");
    %save(['.\Results_' ptid '\res_' ptid '_CARD.mat'],'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
    %save(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_CARD_' int2str(N_iter) '.mat'],'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
    %save(['.\Results_' ptid '\res_' ptid '_CARD_60.mat'],'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
end