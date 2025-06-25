%clc;clear;close all;
ptid='9306087'; 
gridDim = '111';
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=1;
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);

Rx = 4;
Ry = 8;
Rz = 8;
dr = 10;
ddr = 3;
% Rx = [4];%0:2:8;%0:2:8;%0;%0:2:8;
% Ry = 0:2:8;%3;%0:2:8;
% Rz = 0:2:8;%4;%0:2:8;
for Ri = Rx
for Rj = Ry
for Rk = Rz
R = [Ri,Rj,Rk]; 
%load([folder 'MBRT_lattice\LVertices113_' ptid '.mat']);
%fname = strcat('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\LVertices113_9306087_',num2str(R(1)),num2str(R(2)),num2str(R(3)),'.mat');
fname = strcat('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\LVertices\LVertices' ,gridDim, '_',num2str(dr),'_',num2str(ddr),ptid,num2str(R(1)),num2str(R(2)),num2str(R(3)),'.mat');
if exist(fname)
load(fname);
else 
    continue
end
ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy; Dmax<54Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy; Dmax<36Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy; Dmax<54Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy; Dmax<54Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc; Dmax<60Gy


N_iter=30;

Vpeak = setdiff(Vpeak,Vlattice{3}{1});
Vvalley = [Vvalley;Vlattice{3}{1}];

mod_ddr = 5;
mod_Vvalley = Vvalley;
for i = 1:size(Plattice,1)
    x=Plattice(i,1);
    y=Plattice(i,2);
    z=Plattice(i,3);
    id = generate_id(x,y,z,mod_ddr,doseGrid);
    mod_Vvalley = setdiff(mod_Vvalley,id);
end

ctv=cell(1,1);
ctv{1}=ctv1;
n_oar=zeros(N_oar,1);
for i=1:N_oar
    n_oar(i)=numel(oar{i});
end
c=[{Vpeak};{mod_Vvalley};{body};oar;];
%c=[{Vpeak};{mod_Vvalley};{setdiff(Vvalley,mod_Vvalley)};{body};oar;];

N_c=numel(c);
n_c=zeros([N_c 1]);
for i=1:N_c
    n_c(i)=numel(c{i});
end

px=2;nfrac=10;px0=px*nfrac;pvdr = 5;

N_obj = 13;
w_obj = [1; 100; 1; 1; 10; 100; 0.1; 1; 1; 1; 1; 1; 1;];
s_obj = [pvdr*px; pvdr*px; pvdr*px * 1.1; px; px; px * 1.1; 0; px; [15; 10; 10; 10; 12]/nfrac];
n_obj = round([nan; n_c(1) * 0.95; nan; nan; n_c(1) * 0.95; nan; nan; nan; nan; nan; nan; nan; 5 / 0.3^3;]);
c_obj = [1; 1; 1; 2; 2; 2; 3; 3; 4; 5; 6; 7; 8;];
id_obj = cell(N_obj, 1);
type_obj = [0; 2; 3; 0; 2; 3; 0; 3; 3; 3; 3; 3; 1;];

% N_obj = 16;
% w_obj = [1; 100; 1; 1; 1; 100; 1; 1; 1; 0.1; 1; 1; 1; 1; 1; 1;];
% s_obj = [pvdr*px; pvdr*px; pvdr*px * 1.1; px; px; px * 1.1; px; px; px * 1.1; 0; px; [15; 10; 10; 10; 12]/nfrac];
% n_obj = round([nan; n_c(1) * 0.95; nan; nan; n_c(1) * 0.95; nan; nan; n_c(1) * 0.95; nan; nan; nan; nan; nan; nan; nan; 5 / 0.3^3;]);
% c_obj = [1; 1; 1; 2; 2; 2; 3; 3; 3; 4; 4; 5; 6; 7; 8; 9;];
% id_obj = cell(N_obj, 1);
% type_obj = [0; 2; 3; 0; 2; 3; 0; 2; 3; 0; 3; 3; 3; 3; 3; 1;];

% 2. load dij
%id=[90 180]; % (depend on setup)
id=[45 135 225 315];
nN=numel(id);
%id=[0 120];
ctc = [7]; %[3 5 7];
coll_angle = [0 90];


load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(1)) '_ctc' num2str(ctc) 'mm_' num2str(id(1)) '.mat']);
Dij1 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(1)) '_ctc' num2str(ctc) 'mm_' num2str(id(2)) '.mat']);
Dij2 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(2)) '_ctc' num2str(ctc) 'mm_' num2str(id(1)) '.mat']);
Dij3 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(2)) '_ctc' num2str(ctc) 'mm_' num2str(id(2)) '.mat']);
Dij4 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(1)) '_ctc' num2str(ctc) 'mm_' num2str(id(3)) '.mat']);
Dij5 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(1)) '_ctc' num2str(ctc) 'mm_' num2str(id(4)) '.mat']);
Dij6 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(2)) '_ctc' num2str(ctc) 'mm_' num2str(id(3)) '.mat']);
Dij7 = dij.physicalDose{1};
load([folder ptid '\dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(2)) '_ctc' num2str(ctc) 'mm_' num2str(id(4)) '.mat']);
Dij8 = dij.physicalDose{1};
Dij={[Dij1,Dij2,Dij3,Dij4,Dij5,Dij6,Dij7,Dij8]};
%Dij={[Dij1,Dij2,Dij5,Dij6]};
[nY,nX]=size(Dij{1});



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

% % 3. optimization
% load([ptid '\' ptid '_' num2str(numel(id)) '_tvm.mat'])
% load([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'])

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

% wp=struct('Dij',{Dij},'intp',intp,'nX',nX,'nY',nY,'tvx',tvx,'tvz',tvz,...
%     'N_dij',N_dij,'isC',uint32(0));
wp=struct('Dij',{Dij},'nX',nX,'nY',nY,...
    'N_dij',N_dij,'isC',uint32(0));
% wp=struct('Dij',{Dij},'nX',nX,'nY',nY,...
%     'N_dij',N_dij,'isC',uint32(0));
wps=struct('isC',uint32(0));



if mod(method,2)==1
        var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(0),...
            'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
            'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
            'wps',wps,'mu_xs',[],'Wxr',@BX_v1,'Wtxr',@BtX_v1,'wpr',wp,'mu_xr',1);
        tic;x0=admm_mmu1(ip,var_CG);toc;
else
        var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(1),...
            'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
            'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
            'wps',wps,'mu_xs',[],'Wxr',@BX_v3,'Wtxr',@BtX_v3,'wpr',wp,'mu_xr',1);
        tic;x0=admm_mmu4(ip,var_CG);toc;
end

% 4. evaluation
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

%% Get output parameters
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

% % Calculate mean doses in target, body, OAR per angle
% tmp_index = 0;
% mean_doses = zeros(numel(id),N_c);
% for i = 1:numel(id)
%     dd = Dij{1}(:,tmp_index+1:tmp_index+nBeams(i))*double(x(tmp_index+1:tmp_index+nBeams(i)));
%     tmp_index = tmp_index+nBeams(i);
%     for j = 1:N_c
%         mean_doses(i,j) = mean(dd(c{j}));
%     end
% end
% disp(mean_doses);

%% Generate 3D image
% if method>1
d3d=reshape(d,doseGrid.dimensions);
% else
% d=reshape(d,ct.cubeDim);
% end
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
%outp.mean_doses = mean_doses;
%outp.y_bin = reshape(y1,[numel(ctc),numel(id)]);
outp.ObjFnVal = ObjFnVal;

% for i=1:numel(id)
%     if id(i)==0
%         dk=intp_0*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_0 = D_bev;
%         outp.PVDR_0 = [pvdr1,pvdr2];
%     elseif id(i) == 120
%         dk=intp_120*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_120 = D_bev;
%         outp.PVDR_120 = [pvdr1,pvdr2];
%     else
%         dk=intp_240*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_240 = D_bev;
%         outp.PVDR_240 = [pvdr1,pvdr2];
%     end
% end



if mod(method,2)==1
    fname = strcat('output_', ptid, '\res_' ,ptid, '_', num2str(numel(id)), '_', num2str(method), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
    save(fname,"-struct",'outp');
    %save([ptid '\res_' ptid '_' 'tmp.mat'],"-struct",'outp');
    %save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '.mat'],'x','d','d3d','factor', 'mean_doses');
else
    outp.mu_i = mu_i;
    outp.mu_x = mu_x;
    outp.mu_z = mu_z;
    outp.wr = wr;
    save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(numel(ctc)) '_' num2str(method) '.mat'],"-struct",'outp');
end
end
end
end
