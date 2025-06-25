ptid='9306087'; 
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=1;
gridDim = '111';
dr = 10;
ddr = 3;
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);
load('9306087_c.mat')

% 10,4
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-27-19-26.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-29-08-03.mat')

% 10, 6
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-28-07-56.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-10-29-08-02.mat')

% 10,3
%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_2_1_2024-11-04-11-04.mat')
load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-13-02-04.mat')

% ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
% body = cst{1,4}{1};
% N_oar = 5;
% oar = cell(N_oar,1);
% oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy; Dmax<54Gy
% oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy; Dmax<36Gy
% oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy; Dmax<54Gy
% oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy; Dmax<54Gy
% oar(5) = cst{44,4}(1); % Brain V12<5cc; Dmax<60Gy

R = [ 4 8 8 ];
fname = strcat('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\LVertices\LVertices' ,gridDim, '_',num2str(dr),'_',num2str(ddr),'_9306087_',num2str(R(1)),num2str(R(2)),num2str(R(3)),'.mat');
load(fname);

% ctv=cell(1,1);
% ctv{1}=ctv1;
% n_oar=zeros(N_oar,1);
% for i=1:N_oar
%     n_oar(i)=numel(oar{i});
% end
% c=[{Vpeak};{Vvalley};{body};oar;];
Vpeak = setdiff(Vpeak,Vlattice{3}{1});
Vvalley = [Vvalley;Vlattice{3}{1}];
c{1} = Vpeak;
c{2} = Vvalley;

px=2;nfrac=10;px0=px*nfrac;pvdr = 5;

figure; hold on;
for ids = [3 4 8 ]


if ids == 1
d_ADMM = d(c{ids}) / (px*pvdr);
else
d_ADMM = d(c{ids}) / px;
end
N = numel(c{ids});
n = 100;

% Create dose-volume histogram
t = linspace(0, 1.4, n);
dvh_ADMM = zeros(n, 1);
for i = 1:n
    dvh_ADMM(i) = numel(find(d_ADMM >= t(i))) / N;
end

% Plot data
plot(t, dvh_ADMM, 'linewidth', 2);
xlabel('Dose (%)')
ylabel('Volume (%)')
%title('DVH plot for peak and valley doses in the target')
%legend('peak','valley')
 title('DVH plot for OAR')
 legend('body','brainstem','brain')
grid on;


end

hold off;