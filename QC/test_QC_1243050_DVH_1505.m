
ptid='1243050'; 
folder = ['..\' ptid '\'];
addpath('..\utils')
load([folder ptid '.mat'], 'ct', 'cst');
px = 6; % prescription dose
nfrac = 4; % number of fraction


ctv1=cst{11,4}{1};   % ptv 24 Gy, 6 Gy*4
body=cst{1,4}{1};
N_oar=4;
oar=cell(N_oar,1);
oar(1)=cst{2,4}(1); % LargeBowel Dmax<38 Gy, V25<20 cc (D20cc<25 Gy)
oar(2)=cst{3,4}(1); % SmallBowel Dmax<35 Gy, V20<5 cc (D5cc<20 Gy)
oar(3)=cst{12,4}(1); % SpinalCord Dmax<25 Gy
oar(4)=cst{7,4}(1); % L_Kidney 150 cc<12 Gy (D150cc<12 Gy)
oar(5)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};

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
out_folder = "Results_1505_1243050/";
resd = cell(3,1);

% All OAR
fnames = {'Results_1505_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-06-16-09-46.mat' ... %Clinical
    'Results_1505_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-06-16-09-50.mat' ... %GS
    'Results_1505_1243050\res0202_1243050_72_NNZ_3_QC_50_2025-06-14-16-29.mat'}; %QC-BAO
fnames = {'Results_1505_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-06-16-09-46.mat' ... %Clinical
    'Results_1505_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-06-19-11-28.mat' ... %AG
    'Results_1505_1243050\res0202_1243050_72_NNZ_3_QC_50_2025-06-14-16-29.mat'}; %QC-BAO

% fnames = {'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-21-10-04.mat' ... 
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-21-09-53.mat' ...
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-21-09-58.mat'};
% LB+SB
% fnames = {'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-20-13-14.mat' ...
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-20-13-24.mat' ...
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-20-13-21.mat'};
% LB
% fnames = {'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-20-10-21.mat' ...
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-20-11-30.mat' ...
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-20-13-06.mat'};
% No OAR
% fnames = {'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-02-24-10-41.mat' ...
%     'Results_1243050\res0202_1243050_24_NNZ_3_QC_50_2025-03-03-09-26.mat' ...
%     'Results_1243050\res0202_1243050_3_NNZ_3_RND_50_2025-03-03-15-04.mat'};

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
axis([85 115 0 100])
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





