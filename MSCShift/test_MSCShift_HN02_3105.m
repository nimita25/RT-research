
%clc;clear;close all;
ptid = 'HN02';
folder = ['..\' ptid '\'];
addpath('..\utils')
%addpath('G:\Tutorial_030624\code\utils_QC')
load([folder ptid '.mat'], 'ct', 'cst');
px = 8; % prescription dose
nfrac = 5; % number of fraction
mu_min = 5; % MMU threshold
qc = 3; %qc=1: QC, qc=0: mldivide to solve relaxation of MIP, qc = 2: manual choice of shifts, qc = 3: QC with multi MSC per angle

%% Define target and OAR
ctv1 = cst{15,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{9,4}{1};
N_oar = 4;
oar = cell(N_oar,1);
oar(1) = cst{2,4}(1); % R Parotid V50<30Gy
oar(2) = cst{11,4}(1); % OralCavity Dmean<40Gy
oar(3) = cst{17,4}(1); % Oropharynx Dmax<20Gy
oar(4) = cst{16,4}(1); % Larynx Dmax<20Gy


N_Dij = 1;


%clinical choice of 9306087 often used: id=[45,135,225,315]
%% Define  delivery angles
id = [45 135 225 315];

%% Define MSC shifts
if qc <= 1 || qc == 3
    shifts = [0 1 2 -1];

    %% Load influence matrix
    load([folder ptid '_' num2str(id(1)) '_S' num2str(shifts(1)) '.mat'], 'dij');
    m = size(dij.physicalDose{1}, 1);
    N = zeros(numel(id)*numel(shifts), 1);
    N(1) = size(dij.physicalDose{1}, 2);
    if numel(shifts) > 1
        for j = 2:numel(shifts)
            load([folder ptid '_' num2str(id(1)) '_S' num2str(shifts(j)) '.mat'], 'dij');
            N(j) = size(dij.physicalDose{1}, 2);
        end
    end
    
    for i = 2:numel(id)
        for j = 1:numel(shifts)
            load([folder ptid '_' num2str(id(i)) '_S' num2str(shifts(j)) '.mat'], 'dij');
            N((i-1)*numel(shifts)+j) = size(dij.physicalDose{1}, 2);
        end
    end
    
    n_gs = zeros(numel(id)*numel(shifts),1);
    Dij = sparse(m, sum(N));
    n = 0;
    for i = 1:numel(id)
        for j = 1:numel(shifts)
            load([folder ptid '_' num2str(id(i)) '_S' num2str(shifts(j)) '.mat'], 'dij' );
            Dij(:, n + (1:N((i-1)*numel(shifts)+j))) = dij.physicalDose{1};
            n_gs((i-1)*numel(shifts)+j) = N((i-1)*numel(shifts)+j);
            n = n + N((i-1)*numel(shifts)+j);
        end
    end
    [nY, nX] = size(Dij);
    Dij = {Dij};
    
    disp(size(Dij{1}));
    
else
    shifts = [0 0 0 0];

    %% Load influence matrix
    load([folder ptid '_' num2str(id(1)) '_S' num2str(shifts(1)) '.mat'], 'dij');
    m = size(dij.physicalDose{1}, 1);
    N = zeros(numel(id), 1);
    N(1) = size(dij.physicalDose{1}, 2);
    
    for i = 2:numel(id)
        load([folder ptid '_' num2str(id(i)) '_S' num2str(shifts(i)) '.mat'], 'dij');
        N(i) = size(dij.physicalDose{1}, 2);
    end
    
    n_gs = zeros(numel(id),1);
    Dij = sparse(m, sum(N));
    n = 0;
    for i = 1:numel(id)
        load([folder ptid '_' num2str(id(i)) '_S' num2str(shifts(i)) '.mat'], 'dij' );
        Dij(:, n + (1:N(i))) = dij.physicalDose{1};
        n_gs(i) = N(i);
        n = n + N(i);
    end
    [nY, nX] = size(Dij);
    Dij = {Dij};
    
    disp(size(Dij{1}));


    shifts = [0;];
    
end

NNZ = numel(id);
N_iter = 10;
N_iter_admm = 50;



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
N_obj = 9;
type_obj = [0;2;3;0;3;3;3;3;3;];
w_obj = [1;40;40;0.1;1;1;1;1;1;];
s_obj = [px;px;px*1.1;0;px;[30; 40; 20; 20]/nfrac;];
n_obj = round([nan;n_c(1) * 0.95;nan;nan;nan;n_c(3) * 0.50;n_c(4) * 0.50;nan;nan;]);
c_obj = [1;1;1;2;2;3;4;5;6;];


% default values
% N_obj = 9;
% type_obj = [0;2;3;0;3;3;3;3;3;];
% w_obj = [1;40;40;0.1;1;1;1;1;1;];
% s_obj = [px;px;px*1.1;0;px;[30; 40; 20; 20]/nfrac;];
% n_obj = round([nan;n_c(1) * 0.95;nan;nan;nan;n_c(3) * 0.50;n_c(4) * 0.50;nan;nan;]);
% c_obj = [1;1;1;2;2;3;4;5;6;];
% values used in QC-BAO exps
% N_obj = 12;
% type_obj = [0;2;3;0;3;3;3;3;3;0;0;0];
% w_obj = [1;40;40;0.1;1;1;1;1;1;0.1;0.5;0.5];
% s_obj = [px;px;px*1.1;0;px;[30; 40; 20; 20]/nfrac;0;0;0];
% n_obj = round([nan;n_c(1) * 0.95;nan;nan;nan;n_c(3) * 0.50;n_c(4) * 0.50;nan;nan;nan;nan;nan]);
% c_obj = [1;1;1;2;2;3;4;5;6;4;5;6];

load([folder 'dij_' ptid '_doseGrid113.mat']);
for i=1:N_c
tmpCube=zeros(ct.cubeDim);
tmpCube(c{i}) = 1;
% interpolate cube
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
c{i} = VdoseGrid;
end

for i=1:N_c
    n_c(i)=numel(c{i});
end

id_obj = cell(N_obj, 1);
for i=1:N_obj
    id_obj(i,:)=c(c_obj(i));
end
% AmX = @AmX_robust;
% AtmX = @AtmX_robust;
% %AtmX = @AtmX_robust_updated;
% Update_ac = @update_ac_robust;
% Calc_obj_dvh = @calc_obj_dvh_robust;
% % Update_ac = @update_ac;
% % Calc_obj_dvh = @calc_obj_dvh;

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
Ni = 100;
ns = numel(n_gs);
wp = struct('Nx',uint32(nX),'na',uint32(ns),'shifts',shifts,'id',id,'nx',uint32(n_gs),'isC',uint32(0),'w',ones([nX 1],'single'));
ip = struct('Ni',Ni,'N_iter',N_iter,'update_ac',Update_ac,'calc_obj_dvh',Calc_obj_dvh,'isC',uint32(0),'nplot',10000,'var_plot', var_plot);
x0 = ones([nX 1],'single');
maxAtA1 = norm(AtmX(AmX(x0,para0),para0))/norm(x0);


%% Optimization
var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
    'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', [],'wpr',wp);

nnz_x = 10;
ip.mup = maxAtA1*1e-2;%1e-1
disp(maxAtA1)
disp(ip.mup)
ip.mu = 5;%1e-11;%maxAtA2*1e-1;
ip.mu_min = mu_min;
ip.nnz_x = nnz_x;
x_init = zeros([nX 1], 'single');

if qc <= 1 || qc == 3
%tic; [x0,s] = admm_mip_2601(ip,var_CG,NNZ,x_init,qc); toc; %sg := active beam angles

if qc <=1 
tic; [x0,s] = admm_mipMSCShift_3105(ip,var_CG,x_init,qc);  %s := active shifts
time_mip = toc;

n = 0;
tmp_N = zeros(numel(id),1);
y1 = zeros(numel(id)*numel(shifts),1);
for i_id = 1:numel(id)
    [~,tmp_id] = max(s(n+(1:numel(shifts))));
    y1(n+tmp_id) = 1;
    tmp_N(i_id) = N(n+tmp_id);
    n = n+numel(shifts);
end
N = tmp_N;

else
tic; [x0,s] = admm_mipmultiMSCShift_3105(ip,var_CG,x_init,NNZ,qc); %s := active shifts
time_mip = toc;

y1=s;

end

nn = 0;
%var_CG.var_AtA.Dij = [];
tmp_Dij = [];
for ii = 1:numel(id)*numel(shifts)
    if y1(ii) == 1
        tmp_Dij = [tmp_Dij Dij{1}(:,nn+(1:n_gs(ii)))];
    end
    nn = nn+n_gs(ii);
end
Dij = {tmp_Dij};
tmp_Dij = {tmp_Dij};
[nY, nX] = size(tmp_Dij{1});
x_init = zeros([nX 1], 'single');

para.Dij = tmp_Dij;
para.nX = nX;
para.nY = nY;
para0.Dij = tmp_Dij;
para0.nX = nX;
para0.nY = nY;
x0 = ones([nX 1],'single');
maxAtA1 = norm(AtmX(AmX(x0,para0),para0))/norm(x0);
ip.mup = maxAtA1*1e-2;%1e-1
disp(ip.mup)


elseif qc == 2
y1 = ones(numel(id),1);

end

%% Call ADMM method to optimize over active shifts

% mup = 0.01;
% ip.mup = maxAtA * mup;
ip.N_iter = N_iter_admm;

var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
    'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', []);
tic; x1 = admm_mmu_2601(ip, var_CG);
time1 = toc;

x0 = x1;


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
%d3d = reshape(d,ct.cubeDim);

para.Dij = Dij;
[nY, nX] = size(Dij{1});
para.nY = nY;
para.nX = nX;
[tmp,v] = Update_ac(x,para);
[obj,D] = Calc_obj_dvh(x,v);
[mean(sum(obj)) D95 Dmax]
%dlmwrite('GridSearch_SmartInit.csv', [i_mup i_mu mean(sum(obj)) D95 Dmax],'delimiter',',','-append');
obj_total = mean(sum(obj));

[D95, Dmax, CI, Dmean_oar1, Dmax_oar1, Dmean_oar2, Dmax_oar2, Dmean_oar3, Dmax_oar3, Dmean_oar4, Dmax_oar4, Dmean_body] = calcpara_HN02_0202(nfrac, d, px, c);
resd = [D95, Dmax, CI, Dmean_oar1, Dmax_oar1, Dmean_oar2, Dmax_oar2, Dmean_oar3, Dmax_oar3, Dmean_oar4, Dmax_oar4, Dmean_body];

% Save output
clear outp
outp.id = id;
outp.shifts = shifts;
outp.obj_total = obj_total;
outp.para = rmfield(para,'Dij');
outp.para = rmfield(outp.para,'c');
outp.time_admm = time1;
if exist('time_mip') == 1
outp.time_mip = time_mip;
else
outp.time_mip = 0;
end
outp.x = x;
outp.factor = factor;
outp.D95 = D95;
outp.Dmax = Dmax;
outp.d = d;
%outp.d3d = d3d;
outp.n_gs = n_gs;
outp.y1 = y1;
outp.obj = obj;
outp.x0 = x0;
outp.xp = xp;
%outp.s = s;
outp.resd = resd;
outp.mu = ip.mu;
outp.mu_min = mu_min;
if qc == 0
outp.method = 'MIP';
elseif qc == 1
outp.method = 'QC';
elseif qc == 2
outp.method = 'RND';
elseif qc == 3
outp.method = 'multiQC';
end
%fname =strcat('.\Results_', ptid, '\res_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter_admm), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
fname =strcat('.\Results_', ptid, '\res0602_', ptid, '_', int2str(numel(id)), '_shifts_', int2str(numel(shifts)), '_', outp.method, '_', int2str(N_iter_admm), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
save(fname,"-struct","outp");


