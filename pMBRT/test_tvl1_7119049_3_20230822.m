clc;clear;close all;
ptid='7119049';
addpath(genpath('D:\KUMC\matRad-master'));
method=19;
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder=['D:\KUMC\newdata101921\' ptid '\'];
load([folder ptid '_ff.mat'],'ct','cst');

ctv1=cst{9,4}{1};   % ptv60, 2*30
body=cst{1,4}{1};
N_oar=3; % (!)
oar=cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
oar(4)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};
N_iter=30;

ctv=cell(1,1);
ctv{1}=ctv1;
n_oar=zeros(N_oar,1);
for i=1:N_oar
    n_oar(i)=numel(oar{i});
end
c=[ctv;{body};oar;];

N_c=numel(c);
n_c=zeros([N_c 1]);
for i=1:N_c
    n_c(i)=numel(c{i});
end

px=2;nfrac=30;px0=px*nfrac;

N_obj=10;
type_obj=[0;2;3;0;3;
    1;1;1;1;0];
flag_obj=[0;0;0;0;0;
    0;0;0;0;0];
w_obj=[1;10;1;0.1;1;
    1;1;1;1;0.1];
s_obj=[px0;px0;px0*1.1;0;px0;
    [18;12;27;20;0]]*px/px0;
n_obj=round([nan;n_c(1)*0.95;nan;nan;nan;
    n_c(3)*0.5;n_c(3)*0.3;n_c(4)*0.6;n_c(5)*0.5;nan;]);
c_obj=[1;1;1;2;2;
    3;3;4;5;4];

% 2. load dij
id=[0 120 240]; % (depend on setup)
nN=numel(id);
N_dij=1;
tic;
switch method
    case {1}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc0mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc0mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc0mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}*2;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {3,4}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/700;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc7mm_b3_RS10mm_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc9mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}/30;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {5,6}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/600;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc6mm_b3_RS10mm_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc8mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}/40;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {7,8}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/6;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc6mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc8mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}*3;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {9,10}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/6;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc7mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc9mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}*3;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {11,12}
        load([folder 'dij_7119049_collimator_ctv_beam3_0deg+120_ctc3mm_bw3_doseGrid113_0.mat']);
        Dij1=Dij{1};
        load([folder 'dij_7119049_collimator_ctv_beam3_0deg+120_ctc5mm_bw3_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_ctv_beam3_0deg+120_ctc7mm_bw3_doseGrid113_240.mat']);
        Dij3=Dij{1};
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {13,14}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/15;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc5mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1}/4;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc7mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1};
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {15,16}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc3mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}*2;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {17,18}
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc5mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1};
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc5mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1}/2;
        load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc5mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1};
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {19}
        load([folder 'dij_0_3.mat']);
        Dij1=Dij/5;
        load([folder 'dij_120_5.mat']);
        Dij2=Dij;
        load([folder 'dij_240_7.mat']);
        Dij3=Dij;
%         load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc7mm_b3_doseGrid113_0.mat']);
%         Dij1=Dij{1}*2;
%         load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc7mm_b3_doseGrid113_120.mat']);
%         Dij2=Dij{1};
%         load([folder 'dij_7119049_collimator_beam3_0deg+120_ctc7mm_b3_doseGrid113_240.mat']);
%         Dij3=Dij{1}*2;
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
end
toc;
load([folder 'dij_7119049_doseGrid113.mat']);
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
save([ptid '_c.mat'],'c');

% 3. optimization
load([ptid '\' ptid '_' num2str(numel(id)) '_tvm.mat'])
load([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'])

mu_min=0;
eps=mu_min/2;

id_obj=cell(N_obj,N_dij);
for i=1:N_obj
    id_obj(i,:)=c(c_obj(i));
end

AmX=@AmX_v3;
AtmX=@AtmX_v3;
var_plot=struct('n_c',n_c,'dmax',px*1.2);
mu_i=2;mu_x=-4;mu_z=0;wr=0.1; % parameters
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

wp=struct('Dij',{Dij},'intp',intp,'nX',nX,'nY',nY,'tvx',tvx,'tvz',tvz,...
    'N_dij',N_dij,'isC',uint32(0));
wps=struct('isC',uint32(0));

switch method
    case {1,3,5,7,9,11,13,15,17,19}
        var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(0),...
            'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
            'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
            'wps',wps,'mu_xs',[],'Wxr',@BX_v1,'Wtxr',@BtX_v1,'wpr',wp,'mu_xr',1);
        tic;x0=admm_mmu1(ip,var_CG);toc;
    case {2,4,6,8,10,12,14,16,18}
        var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(1),...
            'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
            'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
            'wps',wps,'mu_xs',[],'Wxr',@BX_v3,'Wtxr',@BtX_v3,'wpr',wp,'mu_xr',1);
        tic;x0=admm_mmu4(ip,var_CG);toc;
end

% 4. evaluation
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

y=d(c{1});
y2=sort(y,'descend');
Dmax=y2(1)/px;

% if method>1
d=reshape(d,ct.cubeDim);
% else
% d=reshape(d,ct.cubeDim);
% end
figure;imshow3D(d,[0,px]);

if 1
if mod(method,2)==1
        save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '.mat'],'x','d','factor');
else
        save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '_' num2str(mu_i) '_' num2str(-mu_x) '_' num2str(-mu_z) '_' num2str(wr) '.mat'],'x','d','factor');
end
end
