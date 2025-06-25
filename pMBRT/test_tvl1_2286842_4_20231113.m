%clc;clear;close all;
ptid='2286842'; %%% ctv center (85,87,27) %%%
addpath(genpath('D:\KUMC\matRad-master'));
method=2;
% 1. load ct & define oars & optimization parameter (depend on ptid)
%folder=['D:\KUMC\newdata101921\' ptid '\'];
folder=['C:\Users\nshinde\Desktop\pMBRT\' ptid '\' ptid '\'];
load([folder ptid '.mat'],'ct','cst');

ctv1=cst{22,4}{1};   % ptv69.96, 2.12*33
body=cst{21,4}{1};
N_oar=7; % (!)
oar=cell(N_oar,1);
oar(1)=cst{4,4}(1); % right parotid D50<30
oar(2)=cst{5,4}(1); % left parotid Dmean<26
oar(3)=cst{24,4}(1); % larynx Dmean<32 ??
oar(4)=cst{8,4}(1); % esophagus Dmean<22
oar(5)=cst{3,4}(1); % mandible Dmax<73 Dmean<50
oar(6)=cst{1,4}(1); % oral Dmean<35
oar(7)=cst{2,4}(1); % lip Dmean<15
oar(8)={ctv2ptv_080720(ctv1,20,ct.cubeDim,ct.resolution)};
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

px=2.12;nfrac=33;px0=px*nfrac;

N_obj=14;
type_obj=[0;2;3;0;3;
    1;1;1;1;3;1;1;1;0];
flag_obj=[0;0;0;0;0;
    0;0;0;0;0;0;0;0;0];
w_obj=[1;10;1;0.1;1;
    1;1;1;1;1;1;1;1;0.1];
s_obj=[px0;px0;px0*1.1;0;px0;
    [30;26;32;22;73;50;35;15;0]]*px/px0;
n_obj=round([nan;n_c(1)*0.95;nan;nan;nan;
    n_c(3)*0.5;n_c(4)*0.5;n_c(5)*0.5;n_c(6)*0.5;nan;n_c(7)*0.5;n_c(8)*0.5;n_c(9)*0.5;nan;]);
c_obj=[1;1;1;2;2;
    3;4;5;6;7;7;8;9;5];

% 2. load dij
id=[45 135 225 315]; % (depend on setup)
%id = [45 225];
nN=numel(id);
tic;
switch method
    case {1}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc0mm_b3_doseGrid113_45.mat']);
        Dij1=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc0mm_b3_doseGrid113_135.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc0mm_b3_doseGrid113_225.mat']);
        Dij3=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc0mm_b3_doseGrid113_315.mat']);
        Dij4=Dij{1};
        Dij={[Dij1,Dij2,Dij3,Dij4]};
        [nY,nX]=size(Dij{1});
    case {2,3}
        %load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_45.mat']); 
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_45.mat']); 
        Dij1=Dij{1};
        %load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_135.mat']);
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_135.mat']);
        Dij2=Dij{1};
        %load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_225.mat']);
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_225.mat']);
        Dij3=Dij{1};
        %load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_315.mat']);
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_315.mat']);
        Dij4=Dij{1};
        Dij={[Dij1,Dij2,Dij3,Dij4]};
        [nY,nX]=size(Dij{1});
    case {4,5}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_45.mat']);
        Dij1=Dij{1}/2.5;
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc5mm_b3_doseGrid113_135.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc5mm_b3_doseGrid113_225.mat']);
        Dij3=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_315.mat']);
        Dij4=Dij{1}/2.5;
        Dij={[Dij1,Dij2,Dij3,Dij4]};
        [nY,nX]=size(Dij{1});
    case {6,7}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_45.mat']);
        Dij1=Dij{1}/3;
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc5mm_b3_doseGrid113_135.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc5mm_b3_doseGrid113_225.mat']);
        Dij3=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_315.mat']);
        Dij4=Dij{1}/3;
        Dij={[Dij1,Dij2,Dij3,Dij4]};
        [nY,nX]=size(Dij{1});
    case {8,9}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_45.mat']);
        Dij1=Dij{1}/1.7;
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_135.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_225.mat']);
        Dij3=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_315.mat']);
        Dij4=Dij{1}/1.7;
        Dij={[Dij1,Dij2,Dij3,Dij4]};
        [nY,nX]=size(Dij{1});
    case {10,11}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc7mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}*3; %(! scale 100 for 1cm  range shifter, but for three angles scale 30 maybe better)
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {12,13}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_RS10mm_doseGrid113_0.mat']);
        Dij1=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_RS10mm_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc7mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}/40; %(! scale 100 for 1cm  range shifter, but for three angles scale 30 maybe better)
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {14,15}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/100; %(! scale 100 for 1cm  range shifter)
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_RS10mm_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc6mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}/30; %(! scale 100 for 1cm  range shifter, but for three angles scale 30 maybe better)
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {16,17}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/100; %(! scale 100 for 1cm  range shifter)
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_RS10mm_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc7mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}/30; %(! scale 100 for 1cm  range shifter, but for three angles scale 30 maybe better)
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
    case {18,19}
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_doseGrid113_0.mat']);
        Dij1=Dij{1}/100; %(! scale 100 for 1cm  range shifter)
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc3mm_b3_RS10mm_doseGrid113_120.mat']);
        Dij2=Dij{1};
        load([folder 'dij_' ptid '_collimator_beam4_45deg+90_ctc5mm_b3_doseGrid113_240.mat']);
        Dij3=Dij{1}/30; %(! scale 100 for 1cm  range shifter, but for three angles scale 30 maybe better)
        Dij={[Dij1,Dij2,Dij3]};
        [nY,nX]=size(Dij{1});
end
toc;

nBeams = [size(Dij1,2);size(Dij2,2);size(Dij3,2);size(Dij4,2)];

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
save([ptid '_c.mat'],'c');

% 3. optimization
load([ptid '\' ptid '_' num2str(numel(id)) '_tvm.mat'])
load([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'])

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
mu_i=6;mu_x=-6;mu_z=0;wr=0.1; % parameters
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
%wp=struct('Dij',{Dij},'nX',nX,'nY',nY,...
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
d=reshape(d,doseGrid.dimensions);
% else
% d=reshape(d,ct.cubeDim);
% end
figure;imshow3D(d,[0,px]);


tmp_index = 0;
mean_doses = zeros(numel(id),N_c);
for i = 1:numel(id)
    dd = Dij{1}(:,tmp_index+1:tmp_index+nBeams(i))*double(x(tmp_index+1:tmp_index+nBeams(i)));
    tmp_index = tmp_index+nBeams(i);
    for j = 1:N_c
        mean_doses(i,j) = mean(dd(c{j}));
    end
end
disp(mean_doses);

if 1
if mod(method,2)==1
        save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '.mat'],'x','d','factor', 'mean_doses');
else
        save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '_' num2str(mu_i) '_' num2str(-mu_x) '_' num2str(-mu_z) '_' num2str(wr) '.mat'],'x','d','factor', 'mean_doses');
end
end


