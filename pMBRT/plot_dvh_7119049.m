ptid='7119049'; 
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=1;
gridDim = '113';
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);
load([ptid '_c.mat'],'c');


px=2;%nfrac=10;px0=px*nfrac;pvdr = 5;

figure; hold on;
for ids = [2 ]

load('C:\Users\nshinde\Desktop\pMBRT\7119049\res_7119049_3_3_2_ctc357.mat','d')
d_ADMM = d(c{ids}) / px;

N = numel(c{ids});
n = 100;

% Create dose-volume histogram
t = linspace(0, 1.3, n);
dvh_ADMM = zeros(n, 1);
for i = 1:n
    dvh_ADMM(i) = numel(find(d_ADMM >= t(i))) / N;
end

% Plot data
plot(t, dvh_ADMM, 'linewidth', 2);
hold on 

load('C:\Users\nshinde\Desktop\pMBRT\7119049\res_7119049_3_3_2_ctc335.mat','d')
d_ADMM = d(c{ids}) / px;

N = numel(c{ids});
n = 100;

% Create dose-volume histogram
t = linspace(0, 1.3, n);
dvh_ADMM = zeros(n, 1);
for i = 1:n
    dvh_ADMM(i) = numel(find(d_ADMM >= t(i))) / N;
end

% Plot data
plot(t, dvh_ADMM, 'linewidth', 2);

xlabel('Dose (%)')
ylabel('Volume (%)')
%title('DVH plot for peak and valley doses in the target')
legend('DO','MI-DO')
% title('DVH plot for OAR')
% legend('body','brainstem','brain')
grid on;


end

hold off;