
ptid='7119049'; 
folder = ['..\' ptid '\'];
addpath('..\utils')
load([folder ptid '.mat'], 'ct', 'cst');
px = 2; % prescription dose
nfrac = 30; % number of fraction


ctv1=cst{9,4}{1};   % ptv60, 2*30
body=cst{1,4}{1};
N_oar=3; % (!)
oar=cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
oar(4)={ctv2ptv_080720(ctv1,30,ct.cubeDim,ct.resolution)};


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



nm = [3;4;5];
out_folder = "Results_7119049/";
resd = cell(3,1);

% Updated weights
fnames = {'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-03-22-12-17.mat' ...
    'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-03-31-10-06.mat' ...
    'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-03-22-13-55.mat'};
% fnames = {'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-03-22-12-17.mat' ...
%     'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-03-22-13-41.mat' ...
%     'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-03-22-13-55.mat'};

% Original weights
% fnames = {'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-02-26-11-34.mat' ...
%     'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-02-26-14-17.mat' ...
%     'Results_7119049/res0202_7119049_3_NNZ_3_RND_50_2025-02-26-15-14.mat'};


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





