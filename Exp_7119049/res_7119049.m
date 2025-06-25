clc; clear; close all;

%% Define patient ID and prescription dose
ptid = '7119049';
px = 2;

%% Set up the folder path and load necessary data
folder = ['..\' ptid '\'];
load([folder ptid '.mat'], 'ct', 'cst');

%% Initialize cell array to store region indices
c = cell(6,1);
c{1} = cst{9,4}{1};  % PTV60, 2*30
c{2} = cst{1,4}{1};  % Body
c{3} = cst{2,4}{1};  % Cord
c{4} = cst{7,4}{1};  % Lung (Dmean<18Gy, V20<30%)
c{5} = cst{3,4}{1};  % Esophagus (Dmean<20)
c{6} = cst{4,4}{1};  % Heart (V45<60%)

%% Assign regions to variables
ctv = c{1};
body = c{2};
cord = c{3};
lung = c{4};
esoph = c{5};
heart = c{6};

%% Selection of structure based on treatment planning
choice_structure = 2; % 1 for CTV, 2 for Esophagus
switch choice_structure
    case 1
        ids = 1;
    case 2
        ids = 5;
end


%% Load dosimetric data
load(['.\Results_7119049\' 'res_' ptid '.mat']);

%% Calculate dosimetric parameters
[D95, Dmax, CI, Dmean_lung, V20_lung, Dmean_heart, V30_heart, Dmean_esoph, Dmean_cord, Dmean_body] = calcpara_7119049(d, px, ctv, lung, heart, esoph, cord, body);
resd = [obj_total D95 Dmax CI Dmean_lung V20_lung Dmean_heart V30_heart Dmean_esoph Dmean_cord Dmean_body];
resd

%% Adjust doses for selected structure
d = d(c{ids}) / px;
N = numel(c{ids});

% Create dose-volume histogram
n = 100;
t = linspace(0, 1.3, n);
dvh = arrayfun(@(x) sum(d >= x) / N, t);

% Plotting setup
figure; hold on;
plot(t, dvh, 'b', 'LineWidth', 2);
set(gca, 'FontSize', 18);
grid on;

% Conditional plotting adjustments based on structure choice
if choice_structure == 1
    set(gca, 'XLim', [0.9 1.1], 'YLim', [0 1], ...
        'XTick', 0.9:0.1:1.1, 'XTickLabel', 90:10:110, ...
        'YTick', 0:0.2:1, 'YTickLabel', 0:20:100);
else
    set(gca, 'XLim', [0 0.2], 'YLim', [0 0.5], ...
        'XTick', 0:0.1:0.3, 'XTickLabel', 0:10:30, ...
        'YTick', 0:0.1:0.5, 'YTickLabel', 0:10:50);
end

xlabel('Dose (%)');
ylabel('Volume (%)');


% Save the figure based on the selected structure
switch choice_structure
    case 1
        saveas(gcf, './DVH/CTV_7119049.jpg');
    case 2
        saveas(gcf, './DVH/Esophagus_9306087.jpg');
end