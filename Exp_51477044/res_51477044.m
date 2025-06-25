clc; clear; close all;

%% Define patient ID and prescription dose factor
ptid = '51477044';
px = 1.8;

%% Define file paths and load patient data
folder = ['..\' ptid '\'];
load([folder ptid '.mat'], 'ct', 'cst');

%% Initialize and populate structure index cell array with relevant data
c{1}=cst{10,4}{1};   % ptv45 45.0Gy (25*1.8)
c{2}=cst{18,4}{1};
c{3}=cst{9,4}{1}; % bladder D50: 25Gy; D20: 35Gy.
c{4}=cst{17,4}{1}; % rectum D50: 25Gy; D20: 35Gy; D10: 45Gy.
c{5}=cst{11,4}{1}; % femhead_lt D10: 25Gy.
c{6}=cst{12,4}{1}; % femhead_rt D10: 25Gy.
c{7}=union(c{5}, c{6});
c{8}=cst{14,4}{1}; % penilebulb D50: 25Gy.

%% Assign anatomical regions to variables
ctv = c{1};
body = c{2};
bladder = c{3};
rectum = c{4};
oar1 = c{5}; % femhead_lt
oar2 = c{6}; % femhead_rt
femhead = c{7};
penilebulb = c{8};

%% Selection of structure based on treatment planning
choice_structure = 2; % 1 for CTV, 2 for femhead
switch choice_structure
    case 1
        ids = 1;
    case 2
        ids = 3;
end

%$ Load dosimetric data
load(['.\Results_51477044\' 'res_' ptid '.mat'], 'd', 'obj');

%% Calculate DVH parameters
[D95, Dmax, CI, Dmean_bladder, V90_bladder, Dmean_rectum, V90_rectum, Dmean_femhead, Dmean_penilebulb, Dmean_body] = ...
    calcpara_51477044(d, px, ctv, bladder, rectum, femhead, penilebulb, body);
resd = [D95 sum(obj) Dmax CI Dmean_bladder V90_bladder Dmean_rectum V90_rectum Dmean_femhead Dmean_penilebulb Dmean_body];
resd

%% Adjust doses for selected structure
d = d(c{ids}) / px;
N = numel(c{ids});
n = 100;

% Create dose-volume histogram
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
    set(gca, 'XLim', [0 0.8], 'YLim', [0 0.5], ...
        'XTick', 0:0.2:0.8, 'XTickLabel', 0:20:80, ...
        'YTick', 0:0.1:0.5, 'YTickLabel', 0:10:50);
end

xlabel('Dose (%)');
ylabel('Volume (%)');


% Save the figure based on the selected structure
switch choice_structure
    case 1
        saveas(gcf, './DVH/CTV_9306087.jpg');
    case 2
        saveas(gcf, './DVH/femhead_9306087.jpg');
end