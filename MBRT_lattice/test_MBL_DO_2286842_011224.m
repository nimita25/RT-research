ptid='2286842'; 
gridDim = '111'; % set the dimension of the grid. Default = 3x3x3 mm
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=4; %1: M2, 2: M0, 3: M1, 4: conv
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);


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

N_iter=20;


%% Generate peak and valley
if method <= 3
% Set peak center acc to info from Dij --> max dose in 250 plane at 252,325
lattice_x = [252;262;272;252;252;262;262;272;272;272];
lattice_y = [325;315;305;305;315;305;325;315;325;335];
lattice_z = 250*ones(numel(lattice_x),1);%[79;79;79;79;79;79;79;79;79;79];
dr = 10; %min distance between vertex center
ddr = 1.5; %diameter of the vertex
margin = 5;
else %this defines peak vertices for conv. To do: set values acc to voxels in 228... case
lattice_x = [257;272];
lattice_y = [310;325];
lattice_z = 250*ones(numel(lattice_x),1);
dr = 30;
ddr = 8;
margin = 5;
end

%% Set indices of CTV, peak, valley according to gridDim. Note default gridDim = 3x3x3
ctv1=ctv2ptv_080720(ctv1,margin,ct.cubeDim,ct.resolution);
tmpCube=zeros(ct.cubeDim);
tmpCube(ctv1) = 1;
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
ctv1 = VdoseGrid;
if 0%exist('2286842_lattice.mat','file')
    load('2286842_lattice.mat');
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
    Plattice=[lattice_x ,lattice_y,lattice_z];
    save('2286842_lattice.mat', 'Vpeak','Vvalley','Vlattice','Plattice','dr','ddr');
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


px=2.12;nfrac=33;px0=px*nfrac;pvdr = 5;

N_obj=17;
type_obj=[0;2;3;0;2;3;0;3;1;1;1;1;3;1;1;1;0];
w_obj=[1;10;1;1;10;1;0.1;1;1;1;1;1;1;1;1;1;0.1];
s_obj=[px0*pvdr;px0*pvdr;px0*1.1*pvdr;px0;px0;px0*1.1;0;px0;[30;26;32;22;73;50;35;15;0]]*px/px0;
n_obj=round([nan;n_c(1)*0.95;nan;nan;n_c(2)*0.95;nan;nan;nan;n_c(4)*0.5;n_c(5)*0.5;n_c(6)*0.5;n_c(7)*0.5;nan;n_c(8)*0.5;n_c(9)*0.5;n_c(10)*0.5;nan;]);
c_obj=[1;1;1;2;2;2;3;3;4;5;6;7;8;8;9;10;6];
id_obj=cell(N_obj,1);



%% Load dij
%id=[90 180]; % (depend on setup
% id=[0 45 90 135 180 225 270 315];
id = [45 135 225 315];
%id = [45 135];
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
        if idc <= 2
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
%outp.fname = fname;
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
outp.pvdr1 = pvdr1;
outp.pvdr2 = pvdr2;
outp.ObjFnVal = ObjFnVal;



fname = strcat('output_', ptid, '\res_' ,ptid, '_', num2str(numel(id)), '_', num2str(method), '_', string(datetime('now','Format',"yyyy-MM-dd-HH-mm")), '.mat');
save(fname,"-struct",'outp');

