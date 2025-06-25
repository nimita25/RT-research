%clc; clear; close all;

%% Define patient ID and necessary constants
ptid = '51477044';
px = 1.8; % prescription dose
nfrac = 25; % number of fraction
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst'); % Load contour and structure set data

%% Initialize cell array to store contours of interest
c = cell(7,1);
c{1} = cst{10,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
c{2} = cst{18,4}{1};
c{3} = cst{9,4}{1}; % bladder D50: 25Gy; D20: 35Gy.
c{4} = cst{17,4}{1}; % rectum D50: 25Gy; D20: 35Gy; D10: 45Gy.
c{5} = cst{11,4}{1}; % femhead_lt D10: 25Gy.
c{6} = cst{12,4}{1}; % femhead_rt D10: 25Gy.
c{7} = cst{14,4}{1}; % penilebulb D50: 25Gy.

%% Selection of structure based on treatment planning
choice_structure = 1; % 1 for CTV, 2 for bladder
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
oar4 = c{6}; 
oar5 = c{7}; 

% Plotting setup
%figure; hold on;

%% Load ADMM data
%load(['./Results_' ptid '/res_' ptid '.mat']);
%load('./Results_51477044/res_51477044_ADMM_60.mat');
%[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_51477044(d, px, ctv, oar1, oar2, oar3, oar5, body);
%resd = [obj_total, D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body];
%resd(1,[1,2,3,4,5,end])

%% Adjust doses for selected structure
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
% 
% % Plot data
% %plot(t, dvh_ADMM, 'b', 'linewidth', 2);
% plot(t, dvh_ADMM,'Color','black');
%set(gca, 'FontSize', 18);
%grid on;

num_NE = 25:5:35;
PlotColors = {'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', "#7E2F8E"};
%tmp = ['ADMM'];
tmp = [];
figure; hold on;
xlim([0 1.4])


for NNZ = num_NE
    %% Load CARD data
    %load(['./Results_' ptid '/res_' ptid '_CARD.mat']);
    %load('Results_51477044/res_51477044_CARD_60.mat');
    N_iter = 50;
    fname =strcat('.\Results_', ptid, '\res2308_', ptid, '_NNZ_', int2str(NNZ), '_CARD_', int2str(N_iter), '.mat');
    load(fname);
    %load(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_CARD_' int2str(N_iter) '.mat']);
    %load(['.\Results_9306087\' 'res_' ptid '_CARD.mat']);
    %[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_51477044(d, px, ctv, oar1, oar2, oar3, oar5, body);
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
    hold on;

    %% Load MIP data
    %load(['./Results_' ptid '/res_' ptid '_MIP.mat']);
    %load('Results_51477044/res_51477044_MIP_50_10.mat');
    N_iter = 50;
    fname =strcat('.\Results_', ptid, '\res2908_', ptid, '_NNZ_', int2str(NNZ), '_MIP_', int2str(N_iter), 'new.mat');
    load(fname);
    %load(['.\Results_' ptid '\res_' ptid '_NNZ_' int2str(NNZ) '_MIP_' int2str(N_iter) '.mat']);
    %load('G:\Tutorial_030624\code\Exp_51477044\Results_51477044\res_51477044_NNZ_35_MIP_3_2024-08-25-16-06.mat')
    %[D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_51477044(d, px, ctv, oar1, oar2, oar3, oar5, body);
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
    %plot(t, dvh_CARD);
    plot(t, dvh_CARD,'Color',PlotColors{find(num_NE==NNZ)},'LineStyle','--');

    tmp = [tmp;'CARD'+string(NNZ);'MIP'+string(NNZ)];
end


% Add legend
%legend('ADMM', 'CARD', 'MIP')
%legend('ADMM', 'MIP')
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
        legend(tmp,'location','southwest','FontSize',14)
        saveas(gcf,strcat(out_folder,'MultiNNZ_',ptid,'_',string(num_NE(1)),'_',string(num_NE(length(num_NE))),'.fig'));
        %saveas(gcf, './DVH/CTV_combined.jpg');
    case 2
        legend(tmp,'location','northeast','FontSize',14)
        %saveas(gcf, './DVH/BrainStem_Combined.jpg');
        saveas(gcf,strcat(out_folder,'MultiNNZ_BS_',ptid,'_',string(num_NE(1)),'_',string(num_NE(length(num_NE))),'.fig'));
end
