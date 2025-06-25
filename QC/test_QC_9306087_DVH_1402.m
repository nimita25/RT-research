
ptid='9306087'; 
folder = ['..\' ptid '\'];
addpath('..\utils')
load([folder ptid '.mat'], 'ct', 'cst');
px = 2; % prescription dose
nfrac = 10; % number of fraction


ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc

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



nm = [3;7];
out_folder = "Results_9306087/";
resd = cell(3,1);


% Updated weights + QC
fnames = {'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-16-58.mat' ... %Clinical
    'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-17-01.mat' ... %GS
    'Results_9306087\res0202_9306087_4_NNZ_4_RND_50_2025-03-22-10-50.mat' ... %24 angles
    'Results_9306087\res0202_9306087_4_NNZ_4_RND_50_2025-03-22-10-52.mat'}; %72 angles

% Updated weights
% fnames = {'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-16-58.mat' ... %Clinical
%     'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-16-51.mat' ... %GS
%     'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-16-59.mat' ... %24 angles
%     'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-03-21-17-00.mat'}; %72 angles

% Old original weights
% fnames = {'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-02-19-09-52.mat' ... %Clinical
%     'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-02-20-08-52.mat' ... %GS
%     'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-02-28-10-16.mat' ... %24 angles
%     'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-02-28-15-20.mat'}; %72 angles

%GS: 'Results_9306087/res0202_9306087_4_NNZ_4_RND_50_2025-02-20-08-52.mat' ...

%72: res0202_9306087_4_NNZ_4_RND_50_2025-02-28-12-30.mat'}; %72 angles
%72: res0202_9306087_4_NNZ_4_RND_50_2025-02-19-10-52.mat
%24: res0202_9306087_4_NNZ_4_RND_50_2025-02-21-10-15.mat

% Initialize values for DVH plots
n_s = 100;


for ff = 1:numel(fnames)
    
    % Load output here
    fname = fnames{ff};
    disp('Loading output...')
    load(fname,'d');
    
    
    

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
        plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','yellow');
    elseif ff == 3
        plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','blue');
    else
        plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','green');
    end
end
axis([80 125 0 100])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('Clinical','BAO-GS','QC-BAO-24','QC-BAO-72','Location','northeast','FontSize',14);
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
        plot(t(1:tid)*100,100*D_O{nn}((1:tid),ff),'LineWidth',2,'Color','yellow');
    elseif ff == 3
        plot(t(1:tid)*100,100*D_O{nn}((1:tid),ff),'LineWidth',2,'Color','blue');
    else
        plot(t(1:tid)*100,100*D_O{nn}((1:tid),ff),'LineWidth',2,'Color','green');
    end
end
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('Clinical','BAO-GS','QC-BAO-24','QC-BAO-72','Location','northeast','FontSize',14);
grid on
%title(strcat('BED DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DOAR_0202_',ptid,'.fig'));
hold off
end





