
%clc;clear;close all;
ptid = '1243050';
folder = ['G:\Tutorial_030624\code\' ptid '\'];
addpath('G:\Tutorial_030624\code\utils')
addpath('G:\Tutorial_030624\code\utils_QC')
load([folder ptid '.mat'], 'ct', 'cst');
px = 6; % prescription dose
nfrac = 4; % number of fraction
mu_min = 5; % MMU threshold
qc = 2; %qc=1: QC, qc=0: mldivide to solve relaxation of MIP, qc = 2: manual choice

%% Define target and OAR
ctv1=cst{11,4}{1};   % ptv 24 Gy, 6 Gy*4
body=cst{1,4}{1};
N_oar=4;
oar=cell(N_oar,1);
oar(1)=cst{2,4}(1); % LargeBowel Dmax<38 Gy, V25<20 cc (D20cc<25 Gy)
oar(2)=cst{3,4}(1); % SmallBowel Dmax<35 Gy, V20<5 cc (D5cc<20 Gy)
oar(3)=cst{12,4}(1); % SpinalCord Dmax<25 Gy
oar(4)=cst{7,4}(1); % L_Kidney 150 cc<12 Gy (D150cc<12 Gy)
oar(5)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};


N_Dij = 1;


%clinical choice of 9306087 often used: id=[45,135,225,315]
%% Define  delivery angles
if qc <= 1
    %id = 0:5:359;
    id = 0:15:359;
    %id = 15:30:345;
    % angle = [0];
    angle = [0 30 60];
    num_id = numel(id);
    num_angle = numel(angle);
    id = repmat(id,1,num_angle);
    angle = repelem(angle,num_id);
    NNZ = 3;
    
else
    angle = [0 0 0];
    id = [0 120 240];
    % id = [135 240 315];
    % angle = [0 30 60];
    % id = [285 180 90];
    % % id = [270 180 90];
    % id = [0 150 240];

    % angle selection with QC for updated weights
    angle = [30 30 60];
    id = [180 285 90];

    % AG selection
    % angle = [30 60 60];
    % id = [315 150 300];

    % GS selection
    % angle = [0 60 60];
    % id = [255 300 315];
    
    NNZ = 3;
end

N_iter = 10;
N_iter_admm = 50;


%% Load influence matrix
% if exist(['Dij' num2str(numel(id))],'file')
%     load(['Dij' num2str(numel(id))])
% else
load([folder ptid '_' num2str(id(1)) '_' num2str(angle(1)) '.mat'], 'dij', 'stf');
m = size(dij.physicalDose{1}, 1);
N = zeros(numel(id), 1);
N(1) = size(dij.physicalDose{1}, 2);
for i = 2:numel(id)
    disp(id(i))
    load([folder ptid '_' num2str(id(i)) '_' num2str(angle(i)) '.mat'], 'dij');
    N(i) = size(dij.physicalDose{1}, 2);
end
%N = [2365;2373;2378]; %For optimal set from 24 angles

n_gs = zeros(numel(id),1);
Dij = sparse(m, sum(N));
n = 0;
for i = 1:numel(id)
    disp(id(i))
    load([folder ptid '_' num2str(id(i)) '_' num2str(angle(i)) '.mat'], 'dij' , 'stf');
    Dij(:, n + (1:N(i))) = dij.physicalDose{1};
    n_gs(i) = N(i);
    n = n + N(i);
end
Dij = {Dij};
% end
[nY, nX] = size(Dij{1});


disp(size(Dij{1}));

%save(['Dij' num2str(numel(id))], 'Dij', "n_gs");



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
N_obj=15;
type_obj=[0;2;3;0;3;
    3;1;3;1;3;1;0;0;0;0];
w_obj=[1;10;1;0.05;1;1;1;1;1;1;1;0.5;0.1;0.01;0.01];
s_obj=[px;px;px*1.1;0;px;
    [38;25;35;20;25;12]/nfrac;0;0;0;0];
n_obj=round([nan;n_c(1)*0.95;nan;nan;nan;
    nan;20/0.3^3;nan;5/0.3^3;nan;100/0.3^3;nan;nan;nan;nan]);
c_obj=[1;1;1;2;2;
    3;3;4;4;5;6;3;4;5;6];
id_obj = cell(N_obj, 1);

% Original parameters
% N_obj=11;
% type_obj=[0;2;3;0;3;
%     3;1;3;1;3;1];
% w_obj=[1;10;1;0.1;1;
%     1;1;1;1;1;1];
% s_obj=[px0;px0;px0*1.1;0;px0;
%     [38;25;35;20;25;12]]*px/px0;
% n_obj=round([nan;n_c(1)*0.95;nan;nan;nan;
%     nan;20/0.3^3;nan;5/0.3^3;nan;150/0.3^3]);
% c_obj=[1;1;1;2;2;
%     3;3;4;4;5;6];

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
wp = struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(n_gs),'isC',uint32(0),'w',ones([nX 1],'single'));
ip = struct('Ni',Ni,'N_iter',N_iter,'update_ac',Update_ac,'calc_obj_dvh',Calc_obj_dvh,'isC',uint32(0),'nplot',10000,'var_plot', var_plot);
x0 = ones([nX 1],'single');
maxAtA1 = norm(AtmX(AmX(x0,para0),para0))/norm(x0);


%% Optimization
var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
    'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', [],'wpr',wp);

nnz_x = 10;
ip.mup = maxAtA1*1e-2;%1e-1
ip.mu = 1e2;%1e-11;%maxAtA2*1e-1;
ip.mu_min = mu_min;
ip.nnz_x = nnz_x;
x_init = zeros([nX 1], 'single');

if qc <= 1
tic; [x0,s] = admm_mip_2601(ip,var_CG,NNZ,x_init,qc); toc; %sg := active beam angles

[mk, sg] = maxk(s,NNZ);
else
%sg = [2;4;6;8]; %set sg values to choose clinincal angles
s = ones(NNZ,1);
sg = 1:NNZ;
end

%% Call ADMM method to optimize over active beam angles
nn = 0;
%var_CG.var_AtA.Dij = [];
tmp_Dij = [];
for ii = 1:ns
    if ismember(ii,sg)
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

%% Optimization using ADMM method
% mup = 0.01;
% ip.mup = maxAtA * mup;
ip.N_iter = N_iter_admm;

var_CG = struct('id_x12', uint32(0), 'id_xs', uint32(1), 'id_xr', uint32(0),'JtJ', 'JtJ_id', 'cg_iter', 10,...
    'cg_tol', 1e-5,'AmX', AmX, 'AtmX', AtmX, 'var_AtA', para, 'Wxs', [], 'Wtxs', [], 'wps', [], 'mu_xs', []);
x1 = admm_mmu_2601(ip, var_CG);

x0 = x1;
% % Scale up x1 to original size
% nX = size(Dij{1},2);
% x0 = zeros(nX,1);
% nn = 0;
% nn1 = 0;
% for ii = 1:ns
%     if ismember(ii,sg)
%         x0(nn+(1:n_gs(ii))) = x1(nn1+(1:n_gs(ii)));
%         nn1 = nn1+n_gs(ii);
%     end
%     nn = nn+n_gs(ii);
% end

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
% if length(s) == 72
% v.w_obj(1) = 1;
% v.w_obj(4) = 0.092;
% end
[obj,D] = Calc_obj_dvh(x,v);
[mean(sum(obj)) D95 Dmax]
%dlmwrite('GridSearch_SmartInit.csv', [i_mup i_mu mean(sum(obj)) D95 Dmax],'delimiter',',','-append');
obj_total = mean(sum(obj));

[D95, Dmax, CI, Dmax_oar1, V25_oar1, Dmax_oar2, V20_oar2, Dmax_oar3, V25_oar3, Dmax_oar4, V12_oar4, Dmean_body] = calcpara_1243050_0202(nfrac, d, px, c);
resd = [D95, Dmax, CI, Dmax_oar1, V25_oar1, Dmax_oar2, V20_oar2, Dmax_oar3, V25_oar3, Dmax_oar4, V12_oar4, Dmean_body];

% Save output
clear outp
outp.id = id;
outp.angle = angle;
outp.obj_total = obj_total;
outp.para = para;
outp.x = x;
outp.factor = factor;
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
if qc == 0
outp.method = 'MIP';
elseif qc == 1
outp.method = 'QC';
else
outp.method = 'RND';
end
%fname =strcat('.\Results_', ptid, '\res_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter_admm), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
fname =strcat('.\Results_', ptid, '\res0202_', ptid, '_', int2str(numel(id)), '_NNZ_', int2str(NNZ), '_', outp.method, '_', int2str(N_iter_admm), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
%save(fname,'obj_total', 'x', 'x0', 'xp', 'd', 'obj', 'D95', 'Dmax', 'n_gs');
save(fname,"-struct","outp");


