
ptid='HN02'; 
folder = ['..\' ptid '\'];
addpath('..\utils')
load([folder ptid '.mat'], 'ct', 'cst');
px = 8; % prescription dose
nfrac = 5; % number of fraction


ctv1 = cst{15,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{9,4}{1};
N_oar = 4;
oar = cell(N_oar,1);
oar(1) = cst{2,4}(1); % R Parotid V50<30Gy
oar(2) = cst{11,4}(1); % OralCavity Dmean<40Gy
oar(3) = cst{17,4}(1); % Oropharynx Dmax<20Gy
oar(4) = cst{16,4}(1); % Larynx Dmax<20Gy

ctv = cell(1,1);
ctv{1} = ctv1;
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);
    n_oar(i) = numel(oar{i});
end
c = [ctv; {body}; oar;];
%==============================================
N_c = numel(c);
n_c = zeros([N_c 1]);
for i = 1:N_c
    n_c(i) = numel(c{i});
end

nm = [3;4;5;6];
out_folder = "Results_HN02/";
resd = cell(3,1);

% Dmax wt = 40; Dmean for OP, L
fnames = {'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-11-37.mat' ...%Clincal
    'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-31-09-28.mat' ...%GS
    'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-11-39.mat'};%72 angles

% fnames = {'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-11-37.mat' ...%Clincal
%     'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-12-48.mat' ...%24 angles
%     'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-21-11-39.mat'};%72 angles

% No Dmean OAR
% fnames = {'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-05-11-31.mat' ...%Clincal
%     'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-05-11-30.mat' ...%24 angles
%     'Results_HN02\res0202_HN02_4_NNZ_4_RND_50_2025-03-06-15-08.mat'};%72 angles

%24: res0202_1243050_24_NNZ_3_QC_50_2025-03-03-09-26.mat
%72: res0202_1243050_3_NNZ_3_RND_50_2025-03-03-15-04.mat

% Initialize values for DVH plots
n_s = 100;


for ff = 1:numel(fnames)
    
    % Load output here
    fname = fnames{ff};
    disp('Loading output...')
    load(fname);
    
    
    

    % Generating DVH plot for target
    %DVH_eval = BED2(c{1})/BED_target;
    DVH_eval = d(c{1})/px;
    N = numel(DVH_eval);
    % Create dose-volume histogram for target 
    t = linspace(0, 1.2, n_s);
    for i = 1:n_s
        D_T(i,ff) = numel(find(DVH_eval >= t(i))) / N;
    end

    % Generating DVH plot for OAR
    for nn = 1:numel(nm)
        DVH_eval = d(c{nm(nn)})/px;
        N = numel(DVH_eval);
        % Create dose-volume histogram for target BED
        t = linspace(0, 1.2, n_s);
        for i = 1:n_s
            D_O{nn}(i,ff) = numel(find(DVH_eval >= t(i))) / N;
        end
    end


end



% Generate plots

%Generate DVH plot for target
figure
hold on
t = linspace(0, 1.2, n_s);
tid = 100;
for ff = 1:numel(fnames)
    if ff == 1
    plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','red');
    elseif ff == 2
        plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','blue');
    else
        plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','green');
    end
end
axis([80 125 0 100])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('Clinical','QC-BAO-24','QC-BAO-72','Location','northeast','FontSize',14);
grid on
%title(strcat('BED DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DTumor_0202_',ptid,'.fig'));
hold off

% Generate DVH plots for OAR
for nn = 1:numel(nm)
figure
hold on
t = linspace(0, 1.2, n_s);
for ff = 1:numel(fnames)
    if ff == 1
    plot(t(1:tid)*100,100*D_O{nn}((1:tid),ff),'LineWidth',2,'Color','red');
    elseif ff == 2
        plot(t(1:tid)*100,100*D_O{nn}((1:tid),ff),'LineWidth',2,'Color','blue');
    else
        plot(t(1:tid)*100,100*D_O{nn}((1:tid),ff),'LineWidth',2,'Color','green');
    end
end
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('Clinical','QC-BAO-24','QC-BAO-72','Location','northeast','FontSize',14);
grid on
%title(strcat('BED DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DOAR_0202_',ptid,'.fig'));
hold off
end





