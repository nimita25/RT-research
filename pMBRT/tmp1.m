ptid='7119049';
addpath(genpath('D:\KUMC\matRad-master'));
method=2;
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder=['C:\Users\nshinde\Desktop\pMBRT\2286842\' ptid '\'];
load([folder ptid '.mat'],'ct','cst');

ctv1=cst{9,4}{1};   % ptv60, 2*30
body=cst{1,4}{1};
N_oar=3; % (!)
oar=cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
oar(4)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};
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
%id=[0 120];
ctc = [3 5 7]; %[3 5 7];
nN=numel(id);
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

% 3. optimization
load([ptid '\' ptid '_' num2str(numel(id)) '_tvm.mat'])
load([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'])

mu_min=0;
eps=mu_min/2;

mu_i=2;mu_x=-4;mu_z=0;wr=0.1; % parameters
% para=struct('N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),...
%     'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,...
%     's_obj',s_obj,'n_obj',n_obj,'id_obj',,c_obj,'nN',nN,...
%     'mu_i',mu_i,'mu_x',mu_x,'mu_z',mu_z,'wr',wr, ...
%     'id',id,'ctc',ctc);

wp=struct('intp',intp,'tvx',tvx,'tvz',tvz,...
    'isC',uint32(0));

clear outp
fname = 'C:\Users\nshinde\Desktop\pMBRT\2286842\7119049\res_7119049_3_3_2_ctc375.mat';
outp = load(fname);
d = outp.d;

%% Add the following to _MIP file
wr1=[wr,-wr,-wr];
r_i=mu_i/wr1(1);
r_x=mu_x/wr1(2);
r_z=mu_z/wr1(3);
wp.r_i = r_i;
wp.r_x = r_x;
wp.r_z = r_z;
gd = calc_obj_gd(d,wp);
outp.gdFnVal = gd;
outp.TotalFnVal = outp.gdFnVal+outp.ObjFnVal;
%% Stop additions here

disp(gd); 
disp(outp.TotalFnVal)
save(fname,"-struct",'outp');
