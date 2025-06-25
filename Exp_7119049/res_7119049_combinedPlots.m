%clc; clear; close all;

%% Define patient ID and necessary constants
ptid = '7119049';
px = 2; % prescription dose
nfrac = 30; % number of fraction
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst'); % Load contour and structure set data

%% Initialize cell array to store contours of interest
c = cell(5,1);
c{1} = cst{9,4}{1};   % ptv60, 2*30
c{2} = cst{1,4}{1};
c{3} = cst{7,4}{1};
c{4} = cst{4,4}{1};
c{5} = cst{3,4}{1};

%% Selection of structure based on treatment planning
choice_structure = 2;
switch choice_structure
    case 1
        ids = 1;
    case 2
        ids = 3;
end

%% Assign regions to variables
ctv = c{1};
body = c{2};
oar1 = c{3}; 
oar2 = c{4}; 
oar3 = c{5}; 

% Plotting setup
figure; hold on;

%% Load ADMM data
%load(['./Results_' ptid '/res_' ptid '.mat']);
load('./Results_7119049/res_7119049_ADMM_60.mat');
%[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, V12_oar5, Dmean_body] = calcpara_7119049(d, px, ctv, oar1, oar2, oar3, body);
%resd = [obj_total, D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body];
%resd(1,[1,2,3,4,5,end])

%% Adjust doses for selected structure
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
plot(t, dvh_ADMM, 'b', 'linewidth', 2);
set(gca, 'FontSize', 18);
grid on;

%% Load CARD data
%load(['./Results_' ptid '/res_' ptid '_CARD.mat']);
load('./Results_7119049/res_7119049_CARD_60.mat');
%load(['.\Results_9306087\' 'res_' ptid '_CARD.mat']);
%[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, V12_oar5, Dmean_body] = calcpara_7119049(d, px, ctv, oar1, oar2, oar3, body);
%resd = [obj_total, D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body];
%resd(1,[1,2,3,4,5,end])

%% Adjust doses for selected structure
d_CARD = d(c{ids}) / px;
N = numel(c{ids});
n = 100;

% Create dose-volume histogram
t = linspace(0, 1.3, n);
dvh_CARD = zeros(n, 1);
for i = 1:n
    dvh_CARD(i) = numel(find(d_CARD >= t(i))) / N;
end

% Plot data
plot(t, dvh_CARD, 'g', 'linewidth', 2);

%% Load MIP data
%load(['./Results_' ptid '/res_' ptid '_MIP.mat']);
load('./Results_7119049/res_7119049_MIP_50_10.mat');
%[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, V12_oar5, Dmean_body] = calcpara_7119049(d, px, ctv, oar1, oar2, oar3, body);
%resd = [obj_total, D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body];
%resd(1,[1,2,3,4,5,end])

%% Adjust doses for selected structure
d_CARD = d(c{ids}) / px;
N = numel(c{ids});
n = 100;

% Create dose-volume histogram
t = linspace(0, 1.3, n);
dvh_CARD = zeros(n, 1);
for i = 1:n
    dvh_CARD(i) = numel(find(d_CARD >= t(i))) / N;
end

% Plot data
plot(t, dvh_CARD, 'r', 'linewidth', 2);


% Add legend
legend('ADMM', 'CARD', 'MIP')
%legend('ADMM', 'MIP')




% Conditional plotting based on structure choice
if choice_structure == 1
    set(gca, 'XLim', [0.9 1.1], 'YLim', [0 1], ...
        'XTick', 0.9:0.1:1.1, 'XTickLabel', 90:10:110, ...
        'YTick', 0:0.2:1, 'YTickLabel', 0:20:100);
else
    set(gca, 'XLim', [0 1], 'YLim', [0 0.7], ...
        'XTick', 0:0.2:1, 'XTickLabel', 0:20:100, ...
        'YTick', 0:0.1:0.7, 'YTickLabel', 0:10:70);
end

xlabel('Dose (%)');
ylabel('Volume (%)');

% Save the figure based on the selected structure
switch choice_structure
    case 1
        saveas(gcf, './DVH/CS1_combined.jpg');
    case 2
        saveas(gcf, './DVH/CS2_Combined.jpg');
end
