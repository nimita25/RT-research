clc;clear;close all;
ptid='1243050'; %%% ctv center (85,87,27) %%%
addpath(genpath('D:\KUMC\matRad-master'));
method=2;
% 1. load ct & define oars & optimization parameter (depend on ptid)
%folder=['D:\KUMC\newdata101921\1243050\'];
folder=['C:\Users\nshinde\Desktop\pMBRT\2286842\' ptid '\'];
load([folder ptid '.mat'],'ct','cst');

ctv1=cst{11,4}{1};   % ptv 24 Gy, 6 Gy*4
body=cst{1,4}{1};
N_oar=4;
oar=cell(N_oar,1);
oar(1)=cst{2,4}(1); % LargeBowel Dmax<38 Gy, V25<20 cc (D20cc<25 Gy)
oar(2)=cst{3,4}(1); % SmallBowel Dmax<35 Gy, V20<5 cc (D5cc<20 Gy)
oar(3)=cst{12,4}(1); % SpinalCord Dmax<25 Gy
oar(4)=cst{7,4}(1); % L_Kidney 150 cc<12 Gy (D150cc<12 Gy)
oar(5)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};
N_iter_mip=5;
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

px=6;nfrac=4;px0=px*nfrac;

N_obj=11;
type_obj=[0;2;3;0;3;
    3;1;3;1;3;1];
w_obj=[1;10;1;0.1;1;
    1;1;1;1;1;1];
s_obj=[px0;px0;px0*1.1;0;px0;
    [38;25;35;20;25;12]]*px/px0;
n_obj=round([nan;n_c(1)*0.95;nan;nan;nan;
    nan;20/0.3^3;nan;5/0.3^3;nan;150/0.3^3]);
c_obj=[1;1;1;2;2;
    3;3;4;4;5;6];

% 2. load dij
id=[0 120 240]; % (depend on setup)
%id=[0 120];
ctc = [3 5 7]; %[3 5 7];
nN=numel(id);
all_Dij = cell(nN,numel(ctc));
mean_Dij = zeros(nN,numel(ctc));
tic;
for i_id=1:nN
    for i_ctc=1:numel(ctc)
        load([folder 'dij_' ptid '_collimator_ctc' num2str(ctc(i_ctc)) 'mm_' num2str(id(i_id)) '.mat']);
        all_Dij(i_id,i_ctc) = dij.physicalDose;
        mean_Dij(i_id,i_ctc) = mean(dij.physicalDose{1}(:));
    end
end
min_mean = min(mean_Dij(:));
for i_id=1:nN
    for i_ctc=1:numel(ctc)
        all_Dij{i_id,i_ctc} = (all_Dij{i_id,i_ctc}/mean(all_Dij{i_id,i_ctc}(:)))*min_mean;
        %disp(mean(all_Dij{i_id,i_ctc}(:)))
    end
end

Dij = [];
for i_id=1:nN
    for i_ctc=1:numel(ctc)
        Dij = [Dij all_Dij{i_id,i_ctc}];
    end
end
Dij = {Dij};
[nY,nX]=size(Dij{1});
toc;

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
mu_i=1;mu_x=-2;mu_z=0;wr=0.1; % parameters
para=struct('Dij',{Dij},'nX',nX,'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),...
    'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,'N_dij',N_dij,...
    's_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'c_obj',c_obj,'nN',nN,...
    'mu_i',mu_i,'mu_x',mu_x,'mu_z',mu_z,'wr',wr, ...
    'all_Dij',{all_Dij},'id',id,'ctc',ctc);
ip=struct('N_iter',[],'nplot',10000,'var_plot',var_plot,'isC',uint32(0),'mu_min',mu_min);


x0=ones([nX 1],'single');
maxAtA=norm(AtmX(AmX(x0,para),para))/norm(x0);

mup=0.01;
mu = 1e-5;
ip.mup=maxAtA*mup;
ip.mu = mu;
ip.N_iter_mip=N_iter_mip;
ip.N_iter=N_iter;

% wp=struct('Dij',{Dij},'intp',intp,'nX',nX,'nY',nY,'tvx',tvx,'tvz',tvz,...
%     'N_dij',N_dij,'isC',uint32(0));
wp=struct('Dij',{Dij},'intp',intp,'nX',nX,'nY',nY,'tvx',tvx,'tvz',tvz,...
    'N_dij',N_dij,'isC',uint32(0));
% wp=struct('Dij',{Dij},'nX',nX,'nY',nY,...
%     'N_dij',N_dij,'isC',uint32(0));
wps=struct('isC',uint32(0));

if mod(method,2)==1
    var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(0),...
            'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
            'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
            'wps',wps,'mu_xs',[],'Wxr',@BX_v1,'Wtxr',@BtX_v1,'wpr',wp,'mu_xr',1);
        %tic;[x0,y0]=admm_mmu1_mip(ip,var_CG);toc;
        n = 0;
        Dij = [];
        nBeams = zeros(numel(ctc),1);
        y1 = zeros(numel(id)*numel(ctc),1);
        for i_id = 1:numel(id)
            %[~,tmp_id] = max(y0(n+(1:numel(ctc))));
            if i_id == 1 
                tmp_id = 1;
            elseif i_id == 2
                tmp_id = 2;
            else
                tmp_id = 3;
            end
            y1(n+tmp_id) = 1;
            Dij = [Dij var_CG.var_AtA.all_Dij{i_id,tmp_id}];
            nBeams(i_id) = size(var_CG.var_AtA.all_Dij{i_id,tmp_id},2);
            n = n+numel(ctc);
        end
        disp(y1);
        Dij = {Dij};
        [nY,nX]=size(Dij{1});
        disp([nY,nX]);
        var_CG.var_AtA.Dij = Dij;
        var_CG.var_AtA.nY = nY;
        var_CG.var_AtA.nX = nX;
        tic;x0=admm_mmu1(ip,var_CG);toc;
else
        var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(1),...
            'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
            'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
            'wps',wps,'mu_xs',[],'Wxr',@BX_v3,'Wtxr',@BtX_v3,'wpr',wp,'mu_xr',1);
        %tic;[x0,y0]=admm_mmu4_mip(ip,var_CG);toc;
        n = 0;
        Dij = [];
        nBeams = zeros(numel(ctc),1);
        y1 = zeros(numel(id)*numel(ctc),1);
        for i_id = 1:numel(id)
            %[~,tmp_id] = max(y0(n+(1:numel(ctc))));
            if i_id == 1 
                tmp_id = 3;
            elseif i_id == 2
                tmp_id = 1;
            else
                tmp_id = 3;
            end
            y1(n+tmp_id) = 1;
            Dij = [Dij var_CG.var_AtA.all_Dij{i_id,tmp_id}];
            nBeams(i_id) = size(var_CG.var_AtA.all_Dij{i_id,tmp_id},2);
            n = n+numel(ctc);
        end
        disp(y1);
        Dij = {Dij};
        [nY,nX]=size(Dij{1});
        disp([nY,nX]);
        var_CG.var_AtA.Dij = Dij;
        var_CG.var_AtA.nY = nY;
        var_CG.var_AtA.nX = nX;
        var_CG.wpr.Dij = Dij;
        var_CG.wpr.nY = nY;
        var_CG.wpr.nX = nX;
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

%% Get output parameters
% Calculate objective function value
[obj,D]=calc_obj_dvh(x,var_CG.var_AtA);
ObjFnVal = mean(sum(obj));
wr1=[wr,-wr,-wr];
r_i=mu_i/wr1(1);
r_x=mu_x/wr1(2);
r_z=mu_z/wr1(3);
wp.r_i = r_i;
wp.r_x = r_x;
wp.r_z = r_z;
gd = calc_obj_gd(d,wp);
outp.gdFnVal = gd;

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
    MD(i_N) = mean(y)*100/px;
end

% Calculate mean doses in in target, body, OAR per angle
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

%% Generate 3D image
% if method>1
d3d=reshape(d,doseGrid.dimensions);
% else
% d=reshape(d,ct.cubeDim);
% end
%figure;imshow3D(d3d,[0,px]);
%% Save output
outp.x = x;
outp.d = d;
outp.d3d = d3d;
outp.factor = factor;
outp.Dmax = Dmax;
outp.CI = CI;
outp.MD = MD;
outp.mean_doses = mean_doses;
outp.ctc = ctc;
outp.id = id;
outp.y_bin = reshape(y1,[numel(ctc),numel(id)]);
outp.ObjFnVal = ObjFnVal;

for i=1:numel(id)
    if id(i)==0
        dk=intp_0*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_0 = D_bev;
        outp.PVDR_0 = [pvdr1,pvdr2];
    elseif id(i) == 120
        dk=intp_120*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_120 = D_bev;
        outp.PVDR_120 = [pvdr1,pvdr2];
    else
        dk=intp_240*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_240 = D_bev;
        outp.PVDR_240 = [pvdr1,pvdr2];
    end
end


% for i=1:numel(id)
%     if id(i)==0
%         dk=intp_0*d;
%         D_0 = mean(dk)*100/px;
%         outp.D_0 = D_0;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_0 = PVDR;
%     elseif id(i) == 120
%         dk=intp_120*d;
%         D_120 = mean(dk)*100/px;
%         outp.D_120 = D_120;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_120 = PVDR;
%     else
%         dk=intp_240*d;
%         D_240 = mean(dk)*100/px;
%         outp.D_240 = D_240;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_240 = PVDR;
%     end
% end

if mod(method,2)==1
    save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(numel(ctc)) '_' num2str(method) '.mat'],"-struct",'outp');
    %save([ptid '\res_' ptid '_' 'tmp.mat'],"-struct",'outp');
    %save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '.mat'],'x','d','d3d','factor', 'mean_doses');
else
    outp.mu_i = mu_i;
    outp.mu_x = mu_x;
    outp.mu_z = mu_z;
    outp.wr = wr;
    save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(numel(ctc)) '_' num2str(method) '.mat'],"-struct",'outp');
end
