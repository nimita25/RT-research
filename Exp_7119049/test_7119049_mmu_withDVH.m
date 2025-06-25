clc;
clear;
close all;

%% Define data path
ptid = '7119049';
folder = ['..\' ptid '\'];
load([folder ptid '.mat'], 'ct', 'cst');
px = 2; % prescription dose
nfrac = 30; % number of fraction
mu_min = 5; % MMU threshold
%% Define target and OAR
ctv1 = cst{9,4}{1};   % ptv60, 2*30
body = cst{1,4}{1};
N_oar = 3;
oar = cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20

%% Define  delivery angles
if 1
    id = [0 120 240];
    N_iter = 60;
else
    id = 0:15:345;
    N_iter = 50;
end

%% Load influence matrix
load([folder ptid '_' num2str(id(1)) '.mat'], 'dij');
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
    load([folder ptid '_' num2str(id(i)) '.mat'], 'dij');
    Dij(:, n + (1:N(i))) = dij.physicalDose{1};
    n = n + N(i);
end
[nY, nX] = size(Dij);
Dij = {Dij};

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
%%
% ==============================================================
%   N_obj     - Number of objectives.
%   w_obj     - Objective weight.
%   s_obj     - Prescription dose.
%   n_obj     - Array stores the number of active index for each objective.
%   c_obj     - Strucuture index for each objective.
%   id_obj    - Active index for each DVH objective.
%   type_obj  - DVH type of each objective.
% ==============================================================
N_obj = 12;
w_obj = [1;1;1;0.1;1;1;1;0.005;1;1;0.01;0.01]; 
s_obj = [px;px;px*1.1;0;px;[18;12;0;27;20;0;0]/nfrac];
n_obj = round([nan;n_c(1)*0.95;nan;nan;nan;n_c(3)*0.5;n_c(3)*0.3;nan;n_c(4)*0.6;n_c(5)*0.5;nan;nan;]);
c_obj = [1;1;1;2;2;3;3;3;4;5;4;5];
id_obj = cell(N_obj, 1);
type_obj = [0;2;3;0;3;1;1;0;1;1;0;0];
AmX = @AmX_v3; 
AtmX = @AtmX_v3;
var_plot = struct('n_c', n_c, 'dmax', px * 1.2);
para = struct('Dij', {Dij}, 'nX', nX, 'nY', nY, 'N_c', N_c, 'n_c', n_c, 'c', {c}, 'isC', uint32(0),...
    'N_obj', N_obj, 'type_obj', type_obj, 'w_obj', w_obj, 's_obj', s_obj, 'n_obj', n_obj, 'id_obj', {id_obj}, 'c_obj', c_obj);
ip = struct('N_iter', [], 'nplot', 10000, 'var_plot', var_plot, 'isC', uint32(0));

%% Initialization
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
maxAtA = norm(AtmX(AmX(x0, para0), para0)) / norm(x0); % One steps of power method to estimate the spectral norm of A^TA
ip.mu_min = mu_min;

%% Optimization
mup = 0.01;
ip.mup = maxAtA * mup;
ip.N_iter = N_iter;
var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
    'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', []);

% ADMM
tic;
x0 = admm_mmu_withDVH(ip, var_CG);
toc;

% Project to MMU 
xp = double(x0);
xp(xp < mu_min / 2) = 0;
xp(intersect(find(xp >= mu_min / 2), find(xp < mu_min))) = mu_min;

%% Plan normalization
% Scaling of x such that D95 = px; 
n_ctv95 = round(n_c(1) * 0.95);
if 1
    y0 = Dij{1} * double(xp);
    y = y0(c{1});
    y2 = sort(y, 'descend');
    factor = px / y2(n_ctv95);
    x = xp * factor;
else
    x = xp;
end

% Plan parameters
d = Dij{1} * double(x);
y = d(c{1});
y2 = sort(y, 'descend');
D95 = y2(n_ctv95) / px;
Dmax = y2(1) / px;
d = reshape(d, ct.cubeDim);

[tmp, v] = update_ac(x, para);
[obj, D] = calc_obj_dvh(x, v);
% plotdvh(D, var_plot);

[mean(sum(obj)), D95, Dmax]
obj_total = mean(sum(obj));
% save(['.\Results_' ptid '\res_' ptid '.mat'],'obj_total', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax');
save(['.\Results_' ptid '\res_' ptid '_ADMM_60.mat'],'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax');
