%clc;clear;close all;
ptid='9306087'; 
gridDim = '111';
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=4; %1: pMBRTL-2, 2: pMBRTL-1, 3: pMBRTL-1e, 4: conv
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);



ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy; Dmax<54Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy; Dmax<36Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy; Dmax<54Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy; Dmax<54Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc; Dmax<60Gy


N_iter=20;

%% Generate peak and valley
if method <= 3
% Set peak center acc to info from Dij
%load('peakid_9306087.mat')
% lattice_x = [248;244;246];
% lattice_y = [228;232;276];
% lattice_x = [252;248;246];
% lattice_y = [232;236;276];
% lattice_x = [248;243;246];
% lattice_y = [228;233;276];
% lattice_x = [258;249;246];
% lattice_y = [238;258;276];
% lattice_z = [79;79;79];
% lattice_x = [249;239;239;249;249];
% lattice_y = [258;248;268;248;268];
% lattice_x = [249;249;249;259;259];
% lattice_y = [258;268;248;258;268];
% lattice_x = [249;249;249;259;259;249;259;259;269;269];
% lattice_y = [258;268;248;258;268;238;248;238;248;268];
% lattice_x = [249;249;259;259;249;259;259;269;269];
% latticmee_y = [268;248;258;268;238;248;238;248;268];
% lattice_x = [249;249;249;259;259;259;259;269;269;269;269;254;254];
% lattice_y = [238;248;268;238;248;258;268;238;248;268;278;263;253];
% lattice_x = [249;249;249;259;259;259;259;269;269;269;269];
% lattice_y = [238;248;268;238;248;258;268;238;248;268;278]; %--> 11 peaks:best so far
lattice_x = [249;249;249;259;259;259;259;269;269;269;269;259;259;]; %13 peaks with best result for 8 angles
lattice_y = [238;248;268;238;248;258;268;238;248;268;278;228;278];
lattice_z = 79*ones(numel(lattice_x),1);%[79;79;79;79;79;79;79;79;79;79];
dr = 10;
ddr = 1.5;
margin = 5;
else
lattice_x = [259;259];
lattice_y = [238;268];
lattice_z = 79*ones(numel(lattice_x),1);
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
    load('9306087_lattice.mat');
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
    save('9306087_lattice.mat', 'Vpeak','Vvalley','Vlattice','Plattice','dr','ddr');
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

px=2;nfrac=10;px0=px*nfrac;pvdr = 5;

N_obj = 13;
w_obj = [1; 100; 1; 1; 1; 100; 0.1; 1; 1; 1; 1; 1; 1;];
s_obj = [pvdr*px; pvdr*px; pvdr*px * 1.1; px; px; px * 1.1; 0; px; [15; 10; 10; 10; 12]/nfrac];
n_obj = round([nan; n_c(1) * 0.95; nan; nan; n_c(2) * 0.95; nan; nan; nan; nan; nan; nan; nan; 5 / 0.3^3;]);
c_obj = [1; 1; 1; 2; 2; 2; 3; 3; 4; 5; 6; 7; 8;];
id_obj = cell(N_obj, 1);
type_obj = [0; 2; 3; 0; 2; 3; 0; 3; 3; 3; 3; 3; 1;];



%% Load dij
%id=[90 180]; % (depend on setup)
id=[0 45 90 135 180 225 270 315];
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
            load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
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
           
            load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
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
        if mod(idc,2) == 0
            load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(0) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
        else
            load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(90) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
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
        load([folder ptid '\dij_' ptid '_' gridDim '_'  num2str(id(idc)) '.mat']);
        if isempty(Dij)
            Dij = {[dij.physicalDose{1}]};
        else
            Dij = {[Dij{1},dij.physicalDose{1}]};
        end
    end
    
end

[nY,nX]=size(Dij{1});


%% Set indices 'c' according to gridDim
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);
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
for i_N = 1:N_c
    y = d(c{i_N});
    MD(i_N) = mean(y);
end
dd = [c{1};c{2}];
[D_bev, pvdr1, pvdr2] = calc_bev_para(d(dd),px);



%% Generate 3D image
d3d=reshape(d,doseGrid.dimensions);
%figure;imshow3D(d3d,[0,px]);

%% Save output
clear outp
outp.gridDim = gridDim;
outp.fname = fname;
outp.dr = dr;
outp.ddr = ddr;
outp.ctc = ctc;
outp.id = id;
outp.w_obj = w_obj;
outp.R = R;
outp.x = x;
outp.d = d;
outp.d3d = d3d;
outp.factor = factor;
outp.Dmax = Dmax;
outp.CI = CI;
outp.MD = MD;
outp.pvdr1 = pvdr1;
outp.pvdr2 = pvdr2;
outp.ObjFnVal = ObjFnVal;



fname = strcat('output_', ptid, '\res_' ,ptid, '_', num2str(numel(id)), '_', num2str(method), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
save(fname,"-struct",'outp');

