%clc;clear;close all;
ptid='1243050'; 
gridDim = '111';
addpath('../pMBRT');
addpath('../utils')
method=1; %1: pMBRTL-2, 2: pMBRTL-1, 3: pMBRTL-1e, 4: conv
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
%load([folder ptid '/' ptid '.mat'],'cst','ct');
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);
clear out_w



%% Define target and OAR
ctv1=cst{11,4}{1};   % ptv 24 Gy, 6 Gy*4
body=cst{1,4}{1};
N_oar=4;
oar=cell(N_oar,1);
oar(1)=cst{2,4}(1); % LargeBowel Dmax<38 Gy, V25<20 cc (D20cc<25 Gy)
oar(2)=cst{3,4}(1); % SmallBowel Dmax<35 Gy, V20<5 cc (D5cc<20 Gy)
oar(3)=cst{12,4}(1); % SpinalCord Dmax<25 Gy
oar(4)=cst{7,4}(1); % L_Kidney 150 cc<12 Gy (D150cc<12 Gy)


N_iter=50;
px=1.8;nfrac=28;px0=px*nfrac;pvdr = 5;

%% Generate peak and valley
if method <= 3
% Set peak center acc to info from Dij --> max dose in 180 z plane at ??
% lattice_x = [291;326];
% lattice_y = [376;386];
% lattice_x = [341;324.9;308.7;292.5;328.9;312.7;296.6;332.9;316.8;300.6];%341
% lattice_y = [355;355;355;335;334;334;334;369;369;369];%324
lattice_x = [341;325;309;292;329;313;297;333;317;300]-1;%341
lattice_y = [355;355;355;355;334;334;334;369;369;369];%324
lattice_z = 180*ones(numel(lattice_x),1);
dr = 10;
ddr = 1.5;
margin = 5;
else
lattice_x = [329;317];
lattice_y = [334;369];
lattice_z = 180*ones(numel(lattice_x),1);
dr = 30;
ddr = 10;
margin = 5;
end

%peakid = setdiff(peakid,[20482154]);

%% Set indices of CTV, peak, valley according to gridDim. Note default gridDim = 3x3x3
ctv1=ctv2ptv_080720(ctv1,margin,ct.cubeDim,ct.resolution);
tmpCube=zeros(ct.cubeDim);
tmpCube(ctv1) = 1;
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
ctv1 = VdoseGrid;
if 0%exist('9306087_lattice.mat','file')
    load('1243050_lattice.mat');
else
    %[lattice_x,lattice_y,lattice_z]  = [peakidx,peakidy,peakidz];% ind2sub(doseGrid.cubeDim,peakid);
    N_obj = numel(lattice_x);
    Vlattice=cell(N_obj,1);
    %disp(N_obj)
    for i=1:N_obj
        Vlattice(i)={{generate_id(doseGrid.x(lattice_x(i)),doseGrid.y(lattice_y(i)),doseGrid.z(lattice_z(i)),ddr,doseGrid)}};
    end
    Vpeak=[Vlattice{:,1}];
    Vpeak=unique(vertcat(Vpeak{:}));
    Vvalley=setdiff(ctv1,Vpeak);
    Plattice=[lattice_x,lattice_y,lattice_z];
    save('1243050_lattice.mat', 'Vpeak','Vvalley','Vlattice','Plattice','dr','ddr');
end



%% Set objective weights
ctv=cell(1,1);
ctv{1}=ctv1;
n_oar=zeros(N_oar,1);
for i=1:N_oar
    n_oar(i)=numel(oar{i});
end
c=[{Vpeak};{Vvalley};{body};oar;];
%c=[{Vpeak};{mod_Vvalley};{setdiff(Vvalley,mod_Vvalley)};{body};oar;];

N_c=numel(c);
n_c=zeros([N_c 1]);
for i=1:N_c
    n_c(i)=numel(c{i});
end



N_obj=14;
type_obj=[0;2;3;0;2;3;0;3;
    3;1;3;1;3;1];
w_obj=[1;10;1;1;10;1;0.1;1;
    1;1;1;1;1;1];
s_obj=[px0*pvdr;px0*pvdr;px0*pvdr*1.1;px0;px0;px0*1.1;0;px0;
    [38;25;35;20;25;12]]*px/px0;
n_obj=round([nan;n_c(1)*0.95;nan;nan;n_c(1)*0.95;nan;nan;nan;
    nan;20/0.3^3;nan;5/0.3^3;nan;150/0.3^3]);
c_obj=[1;1;1;2;2;2;3;3;
    4;4;5;5;6;7];



%% Load dij
%id=[90 180]; % (depend on setup)
id=[0 120 240];
nN=numel(id);
ctc = [7]; %[3 5 7];


if method == 1

    coll_angle = [0 90];
    Dij = {};
    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)
            % if mod(idc,2) == 1
            %     load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(0) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
            % else
            %     load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(90) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
            % end
            load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
            if isempty(Dij)
                Dij = {[dij.physicalDose{1}]};
            else
                Dij = {[Dij{1},dij.physicalDose{1}]};
            end
        end
    end

elseif method == 2
    
    coll_angle = [0];
    Dij = {};
    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)
           
            load([folder  'dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
            if isempty(Dij)
                Dij = {[dij.physicalDose{1}]};
            else
                Dij = {[Dij{1},dij.physicalDose{1}]};
            end
        end
    end
elseif method == 3
    
    Dij = {};
    for idc = 1:numel(id)
        if mod(idc,2) == 1
            load([folder  'dij_' ptid '_' gridDim '_collimator' num2str(0) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
        else
            load([folder  'dij_' ptid '_' gridDim '_collimator' num2str(90) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
        end
        if isempty(Dij)
            Dij = {[dij.physicalDose{1}]};
        else
            Dij = {[Dij{1},dij.physicalDose{1}]};
        end
    end
else
    Dij = {};
    for idc = 1:numel(id)
        load([folder  'dij_' ptid '_' gridDim '_'  num2str(id(idc)) '.mat']);
        if isempty(Dij)
            Dij = {[dij.physicalDose{1}]};
        else
            Dij = {[Dij{1},dij.physicalDose{1}]};
        end
    end
    
end

[nY,nX]=size(Dij{1});


%% Set indices 'c' according to gridDim
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);
for i=3:N_c
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
save([ptid '_c.mat'],'c');

w1 = [10;10;10;10;1;1;1;1];
w2 = [1;10;1;10;1;10;1;10];
w3 = [1;1;10;10;1;1;10;10;];


for w1i = 1:numel(w1)

px=1.8;nfrac=28;px0=px*nfrac;pvdr = 5;
w_obj=[1;w1(w1i);1;1;10;1;0.1;1;
    w2(w1i);1;1;1;1;w3(w1i)];

%% Set hyperparameters
mu_min=0;
eps=mu_min/2;

N_dij=1;
id_obj=cell(N_obj,N_dij);
for i=1:N_obj
    id_obj(i,:)=c(c_obj(i));
end

AmX=@AmX_v3;
AtmX=@AtmX_v3;
var_plot=struct('n_c',n_c,'dmax',px*1.2);
mu_i=1;mu_x=-2;mu_z=0;wr=0.1; % parameters
para=struct('Dij',{Dij},'nX',nX,'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),...
    'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,'N_dij',N_dij,...
    's_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'c_obj',c_obj,'nN',nN,...
    'mu_i',mu_i,'mu_x',mu_x,'mu_z',mu_z,'wr',wr);
ip=struct('N_iter',[],'nplot',10000,'var_plot',var_plot,'isC',uint32(0),'mu_min',mu_min);

x0=ones([nX 1],'single');
maxAtA=norm(AtmX(AmX(x0,para),para))/norm(x0);

mup=0.01;
ip.mup=maxAtA*mup;
ip.N_iter=N_iter;

wp=struct('Dij',{Dij},'nX',nX,'nY',nY,...
    'N_dij',N_dij,'isC',uint32(0));
wps=struct('isC',uint32(0));


%% Optimization
var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(0),...
    'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
    'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
    'wps',wps,'mu_xs',[],'Wxr',@BX_v1,'Wtxr',@BtX_v1,'wpr',wp,'mu_xr',1);
tic;x0=admm_mmu1(ip,var_CG);toc;


%% Calculate output metrics
px = px*pvdr;
xp=double(x0);
xp(xp<eps)=0;
xp(intersect(find(xp>=eps),find(xp<mu_min)))=mu_min;

if 1
n_ctv98=round(n_c(1)*0.95);
y0=Dij{1}*double(xp);
y=y0(c{1});
y2=sort(y,'descend');
factor=px/y2(n_ctv98);
x=xp*factor;
else
x=xp;
end

d=Dij{1}*double(x);

% Get output parameters
% Calculate objective function value
[obj,D]=calc_obj_dvh(x,var_CG.var_AtA);
ObjFnVal = mean(sum(obj));

%Dmax calculation
y=d(c{1});
y2=sort(y,'descend');
Dmax=y2(1)/px; 

%CI calculation
V100 = sum(y2>=px);
V = numel(c{1});
V100_all = sum(d>=px);
CI = (V100*V100)/(V*V100_all); 

%Mean dose in target, body, OAR calculation
MD = zeros(N_c,1);
maxD = zeros(N_c,1);
for i_N = 1:N_c
    y = d(c{i_N});
    MD(i_N) = mean(y);
    maxD(i_N) = max(y);
end
dd = [c{1};c{2}];
[D_bev, pvdr1, pvdr2] = calc_bev_para(d(dd),px);
pvdr3 = mean(d(c{1}))/mean(d(c{2}));



%% Generate 3D image
d3d=reshape(d,doseGrid.dimensions);
%figure;imshow3D(d3d,[0,px]);

%% Save output
clear outp
outp.gridDim = gridDim;
%outp.fname = fname;
% outp.para = rmfield(para,'Dij');
% outp.para = rmfield(outp.para,'c');
outp.lattice_x = lattice_x;
outp.lattice_y = lattice_y;
outp.obj = obj;
outp.dr = dr;
outp.ddr = ddr;
outp.ctc = ctc;
outp.id = id;
outp.w_obj = w_obj;
%outp.R = R;
outp.x = x;
outp.d = d;
outp.d3d = d3d;
outp.factor = factor;
outp.Dmax = Dmax;
outp.CI = CI;
outp.MD = MD;
outp.maxD = maxD;
outp.pvdr1 = pvdr1;
outp.pvdr2 = pvdr2;
outp.pvdr3 = pvdr3;
outp.ObjFnVal = ObjFnVal;

if ~exist('out_w')
    out_w= cell(numel(w1),1);
end
out_w{w1i} = outp;
%counter = counter + 1;
end

fname = strcat('output_', ptid, '/resPareto_' ,ptid, '_', num2str(numel(id)), '_', num2str(method), '.mat');
save(fname,"out_w",'-v7.3');

