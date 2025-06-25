clc;clear;close all;
Initilization
addpath(genpath('D:\luoying\FLASH\matRad-master'));
%%
id = [225 315]; %angle index 
ptid = 'HN02';
px = 40; % prescription dose
N_frac = [5,5,5,5]; % fraction for each field
flag_dro = 1;  
nfrac = 5; % fraction
mup = 0.01; % ρ

%%
N_dij = 9;
for i = 1:N_dij
    load(['D:\luoying\FLASH\IMPT\code\' ptid '\' ptid '_v' num2str(i) '.mat'])
    Dij{i} = Dij0;
end
%%
[nY,nX] = size(Dij{1});
mu_min = 5;
folder = [pwd '\' ptid '\'];
load(['D:\luoying\FLASH\re-irradiation cases\mat\HN02.mat'], 'ct', 'cst');
ctv1 = cst{15,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
% enl=5;
% ptv=ctv2ptv_080720(ctv1,enl,ct.cubeDim,ct.resolution);
% ctv1 = ptv;
body = cst{9,4}{1};
N_oar = 4;
oar = cell(N_oar,1);
oar(1) = cst{2,4}(1); % R Parotid V50<30Gy
oar(2) = cst{11,4}(1); % OralCavity Dmean<40Gy
oar(3) = cst{17,4}(1); % Oropharynx Dmax<20Gy
oar(4) = cst{16,4}(1); % Larynx Dmax<20Gy

for i=1:N_oar
    oar{i} = setdiff(oar{i},ctv1);
end
dx = 3;
shift = 10;
r = ceil(shift/dx);
roi0 = ctv2ptv(ctv1,r,ct.cubeDim);% CTV extension 1cm
roi = setdiff(roi0,ctv1);% CTV1cm
ctv = cell(1,1);
ctv{1} = ctv1;
n_oar = zeros(N_oar,1);% OARs set
for i = 1:N_oar
    n_oar(i) = numel(oar{i});
end
c = [ctv;{body};oar;roi];
N_c = numel(c);
n_c = zeros([N_c 1]);
for i = 1:N_c
    n_c(i) = numel(c{i});
end

N_obj = 9;
type_obj = [0;2;3;0;3;3;3;3;3];
w_obj = [1;40;1;0.1;1;1;1;1;1];
s_obj = [px;px;px*1.1;0;px;[30; 40; 20; 20]];
n_obj = round([nan;n_c(1) * 0.95;nan;nan;nan;n_c(3) * 0.50;n_c(3) * 0.50;nan;nan]);
c_obj = [1;1;1;2;2;3;4;5;6];

id_obj = cell(N_obj, N_dij);
for j = 1:N_dij
    for i = 1:N_obj
        id_obj(i,j) = c(c_obj(i));
    end
end

% %%
% AmX = @AmX_v3;
% AtmX = @AtmX_v3;
AmX = @AmX_Robust;
AtmX = @AtmX_Robust;
plotdvh = @plotdvh_blend1;
var_plot = struct('n_c',n_c,'dmax',px*1.2);
para = struct('Dij',{Dij},'nX',nX,'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),'N_obj',N_obj,'type_obj',type_obj,...
    'w_obj',w_obj,'s_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'N_angle',N_angle,'id_angle',{id_angle},...
    'c_obj',c_obj,'N_frac',N_frac,'N_dij',N_dij,'nfrac',nfrac);

% [m,n] = size(Droi{1});

wp = struct('mu_min',mu_min,'isC',uint32(0),...
   'N_ray',{Nray},'N_angle',N_angle,'id_angle',{id_angle},...
    'N_dij',N_dij,'nfrac',nfrac);
wps = struct('isC',uint32(0));
para.wpr = wp;
N_iter = 30;
var_CG = struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(flag_dro),'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,'AmX',AmX,'AtmX',...
    AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,'wpr',wp,'wps',wps,'mu_xs',[],'mu_xr',[],'nfrac',nfrac);
ip = struct('N_iter',N_iter,'plotdvh',plotdvh,'px',px,'nX',nX,'var_plot',var_plot,'isC',uint32(0),'nplot',N_iter,'mup',[],'nfrac',nfrac);

%% find ρ
x0 = ones([nX 1],'single');
N_obj = 2;
type_obj = [0;0;];
w_obj = [1;0.1;];
s_obj = [px;0;]; % (!)
c_obj = [1;2;]; % (!)
n_obj = round([nan;nan;]);
id_obj = cell(N_obj,N_dij);
for j = 1:N_dij
    for i = 1:N_obj
        id_obj(i,j) = c(c_obj(i));
    end
end
id_obj(1) = c(1); % (!)
id_obj(2) = c(2);
para0 = struct('Dij',{Dij},'nX',nX,'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,...
    's_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'N_angle',N_angle,'id_angle',{id_angle},'c_obj',c_obj,'N_frac',N_frac,'N_dij',N_dij,'nfrac',nfrac);
maxAtA1 = norm(AtmX(AmX(x0,para0),para0))/norm(x0);
ip.mup = maxAtA1*10e-2;
ip.mu_min = mu_min;

%%optimization 
tic;x0 = admm_flash_121823_d_pdr_Robust(ip,var_CG);toc;% tic:start a stopwatch timer; toc:stop the timer

plotdvh = @plotdvh_blend2;
xp = double(x0);
xp(xp < mu_min/2) = 0;                         
xp(intersect(find(xp >= mu_min/2),find(xp < mu_min))) = mu_min;

% Step 2: plan evaluation
N_dij = 9;
px = 40;
% tmp1 = cell(10,9);
% obj1 = cell(12,9); 
%D1 = cell(8,9);
for i = 1:N_dij
    y0{i}= Dij{i}*xp*nfrac;
    Dij_s=Dij{i};
    y{i} = y0{i}(c{1});
    y2{i}  = sort(y{i} ,'descend');
    n_ctv95 = round(n_c(1)*0.95);
    factor_proton{i}  = px/y2{i}(n_ctv95);
    x{i} = xp*factor_proton{i};
    x_s = x{i};
    d{i} = Dij{i}*xp*nfrac;
    d_s = d{i};
    v = para;
    v.Dij = Dij{i};
    [tmp0,v] = update_ac(x_s,v);
    tmp1(:, i) = tmp0;
    [obj0,D0] = calc_obj_dvh(x_s,v);
    %obj1(:, i)= obj0;
    D1(:, i)= D0;
    save(['result/' num2str(ptid) '_result' num2str(i) '.mat'],'Dij0','x_s','xp','x0','d_s','tmp0','obj0','D0');
end
save(['result/' num2str(ptid) '_dose' '.mat'],'D1');
plotdvh(D1,var_plot);


load('D:\luoying\FLASH\IMPT\code\result\HN02_result1.mat');
d2=d_s(ctv1);
Dmax=max(d2)/px;
Dmean=(sum(d2)/numel(ctv1))/px;
CI=sum(d2>=px)^2/(sum(d_s>=px)*numel(ctv1));
