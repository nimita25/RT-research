%clc; clear; close all;

%% Define patient ID and necessary constants
ptid = '9306087';
px = 2;
nfrac = 10;
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst'); % Load contour and structure set data

%% Initialize cell array to store contours of interest
c = cell(7,1);
c{1} = cst{38,4}{1};   % CTV: 20Gy x 1
c{2} = cst{1,4}{1};    % Body
c{3} = cst{37,4}{1};   % Brainstem: Dmax < 15Gy
c{4} = cst{11,4}{1};   % Optic Chiasm: Dmax < 10Gy
c{5} = cst{12,4}{1};   % Optic Nerve Right: Dmax < 10Gy
c{6} = cst{13,4}{1};   % Optic Nerve Left: Dmax < 10Gy
c{7} = cst{44,4}{1};   % Brain: V12 < 5cc

%% Selection of structure based on treatment planning
choice_structure = 2; % 1 for CTV, 2 for Brainstem
switch choice_structure
    case 1
        ids = 1;
    case 2
        ids = 3;
end

%% Assign regions to variables
ctv = c{1};
body = c{2};
oar1 = c{3}; % Brainstem
oar2 = c{4}; % Optic Chiasm
oar3 = c{5}; % Optic Nerve Right
oar4 = c{6}; % Optic Nerve Left
oar5 = c{7}; % Brain

% Plotting setup
figure; hold on;

%% Load ADMM data
%load(['.\Results_9306087\' 'res_' ptid '.mat']);
%load('./Results_9306087/res_9306087_ADMM_60.mat');
%[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_9306087(nfrac, d, px, ctv, oar1, oar2, oar3, oar4, oar5, body);
%resd = [obj_total, D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body];
%resd(1,[1,2,3,4,5,end])

% %% Adjust doses for selected structure
% d_ADMM = d(c{ids}) / px;
% N = numel(c{ids});
% n = 100;
% 
% % Create dose-volume histogram
% t = linspace(0, 1.3, n);
% dvh_ADMM = zeros(n, 1);
% for i = 1:n
%     dvh_ADMM(i) = numel(find(d_ADMM >= t(i))) / N;
% end

% Plot data
%plot(t, dvh_ADMM, 'b', 'linewidth', 2);
% plot(t, dvh_ADMM,'Color','black');
%set(gca, 'FontSize', 18);
grid on;

num_NE = 40:5:50;
PlotColors = {'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', "#7E2F8E"};
%tmp = ['ADMM'];
tmp = [];


for NNZ = num_NE
    %% Load CARD data
    %load(['.\Results_9306087\' 'res_' ptid '_CARD.mat']);
    N_iter = 50;
    fname =strcat('.\Results_', ptid, '\res2308_', ptid, '_NNZ_', int2str(NNZ), '_CARD_', int2str(N_iter), '.mat');
    load(fname);
    %load(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_CARD_' int2str(N_iter) '.mat']);
    %load('./Results_9306087/res_9306087_CARD_60.mat');
    %[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_9306087(nfrac, d, px, ctv, oar1, oar2, oar3, oar4, oar5, body);
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
    %plot(t, dvh_CARD, 'g', 'linewidth', 2);
    plot(t, dvh_CARD,'Color',PlotColors{find(num_NE==NNZ)});
    
    %% Load MIP data
    %load(['.\Results_9306087\' 'res_' ptid '_MIP.mat']);
    N_iter = 50;
    fname =strcat('.\Results_', ptid, '\res2908_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter), 'new.mat');
    load(fname);
    %load(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_MIP_' int2str(N_iter) '.mat']);
    %load('./Results_9306087/res_9306087_MIP_50_10.mat');
    %[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_9306087(nfrac, d, px, ctv, oar1, oar2, oar3, oar4, oar5, body);
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
    %plot(t, dvh_CARD, 'r', 'linewidth', 2);
    plot(t, dvh_CARD,'Color',PlotColors{find(num_NE==NNZ)},'LineStyle','--');

    tmp = [tmp;'CARD'+string(NNZ);'MIP'+string(NNZ)];
end

% Add legend
%legend('ADMM', 'CARD', 'MIP')
disp(tmp)





% Conditional plotting based on structure choice
% if choice_structure == 1
%     set(gca, 'XLim', [0.9 1.1], 'YLim', [0 1], ...
%         'XTick', 0.9:0.1:1.1, 'XTickLabel', 90:10:110, ...
%         'YTick', 0:0.2:1, 'YTickLabel', 0:20:100);
% else
%     set(gca, 'XLim', [0 1], 'YLim', [0 0.7], ...
%         'XTick', 0:0.2:1, 'XTickLabel', 0:20:100, ...
%         'YTick', 0:0.1:0.7, 'YTickLabel', 0:10:70);
% end

xlabel('Dose (%)','FontSize',14);
ylabel('Volume (%)','FontSize',14);

out_folder = "output/";
% Save the figure based on the selected structure
switch choice_structure
    case 1
        %saveas(gcf, './DVH_9306087/CTV_9306087_Combined.jpg');
        legend(tmp,'location','southwest','FontSize',14)
        saveas(gcf,strcat(out_folder,'MultiNNZ_',ptid,'_',string(num_NE(1)),'_',string(num_NE(length(num_NE))),'.fig'));
    case 2
        %saveas(gcf, './DVH_9306087/BrainStem_9306087_Combined.jpg');
        legend(tmp,'location','northeast','FontSize',14)
        saveas(gcf,strcat(out_folder,'MultiNNZ_BS_',ptid,'_',string(num_NE(1)),'_',string(num_NE(length(num_NE))),'.fig'));
end
