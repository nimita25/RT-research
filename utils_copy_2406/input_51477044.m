% ptid='51477044'; % (!)

load([folder ptid '.mat'],'ct','cst');

ctv1=cst{10,4}{1};   % ptv45 45.0Gy (25*1.8)

body=cst{18,4}{1};
N_oar=5; % (!)
oar=cell(N_oar,1);
oar(1)=cst{9,4}(1); % bladder D50: 25Gy; D20: 35Gy.
oar(2)=cst{17,4}(1); % rectum D50: 25Gy; D20: 35Gy; D10: 45Gy.
oar(3)=cst{11,4}(1); % femhead_lt D10: 25Gy.
oar(4)=cst{12,4}(1); % femhead_rt D10: 25Gy.
oar(5)=cst{14,4}(1); % penilebulb D50: 25Gy.

if flag_robust==0
N_dij=1;
% load([folder ptid '_r8_bw5_lss6_0_0_0.mat'],'dij','stf');
load([folder ptid '_r8_bw5_lss3_0_0_0.mat'],'dij','stf');
% load([folder ptid '_proton_0_0_0.mat'],'dij','stf');
Dij=dij.physicalDose(1);
else
N_dij=9;
Dij=cell(N_dij,1);
filename=cell(N_dij,1);
filname{1}='_r8_bw5_lss3_0_0_0.mat';
filname{2}='_r8_bw5_lss3_m5_0_0.mat';
filname{3}='_r8_bw5_lss3_5_0_0.mat';
filname{4}='_r8_bw5_lss3_0_m5_0.mat';
filname{5}='_r8_bw5_lss3_0_5_0.mat';
filname{6}='_r8_bw5_lss3_0_0_m5.mat';
filname{7}='_r8_bw5_lss3_0_0_5.mat';
filname{8}='_r8_bw5_lss3_range_m35.mat';
filname{9}='_r8_bw5_lss3_range_35.mat';
for i=1:N_dij
    load([folder ptid filname{i}],'dij');
    Dij{i}=dij.physicalDose{1};    
end
load([folder ptid filname{1}],'stf');
end

ctv=cell(1,1);
ctv{1}=ctv1;

c=[ctv;{body};oar];
N_c=numel(c);
n_c=zeros([N_c 1]);
for i=1:N_c
    n_c(i)=numel(c{i});
end

nY=size(Dij{1},1);

N_obj=13;
type_obj=[0;2;3;0;3;
    1;1;1;1;1;1;1;1;];
w_obj=[1;100;1;0.1;1;
    1;1;1;1;1;1;1;1;];
s_obj=[45;45;45*1.05;0;45;
    [25;35;25;35;45;25;25;25;]]/25;
n_obj=round([nan;n_c(1)*0.98;nan;nan;nan;
    n_c(3)*0.5;n_c(3)*0.2;
    n_c(4)*0.5;n_c(4)*0.2;n_c(4)*0.1;
    n_c(5)*0.1;n_c(6)*0.1;n_c(7)*0.5;]);
c_obj=[1;1;1;2;2;
    3;3;4;4;4;5;6;7;];
id_obj=cell(N_obj,N_dij);
w_obj=w_obj*ones([1 N_dij]);
para=struct('Dij',[],'N_dij',N_dij,'nX',[],'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),...
    'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,...
    's_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'c_obj',c_obj);

N_obj=2;
type_obj=[0;0;];
w_obj=[1;0.1;];%
s_obj=[45;0;]/25;
n_obj=round([nan;nan;]);
c_obj=[1;2;];
id_obj=cell(N_obj,N_dij);
for i=1:N_obj
id_obj(i,:)=c(i);
end
w_obj=w_obj*ones([1 N_dij]);
para0=struct('Dij',[],'N_dij',N_dij,'nX',[],'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),...
    'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,...
    's_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'c_obj',c_obj);

plotdvh=@plotdvh_prostate_robust;


var_plot=struct('n_c',n_c,'N_dij',N_dij);
