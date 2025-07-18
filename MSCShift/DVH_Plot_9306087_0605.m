
ptid='9306087'; 
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
gridDim = '113';
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['..\' ptid '\'];
load([folder ptid '.mat'], 'ct', 'cst');
load([folder 'dij_' ptid '_doseGrid113.mat']);

%% Define target and OAR
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

for i=1:N_c
tmpCube=zeros(ct.cubeDim);
tmpCube(c{i}) = 1;
% interpolate cube
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
c{i} = VdoseGrid;
end

px = 2; % prescription dose
nfrac = 10; % number of fraction
mu_min = 5; % MMU threshold



clear D_T; clear D_O;

nm = [3;7];
out_folder = "output/";
resd = cell(3,1);
fnames = {'Results_9306087\res0602_9306087_4_shifts_1_RND_50_2025-06-05-17-56.mat' ...
    'Results_9306087\res0602_9306087_4_shifts_4_QC_50_2025-06-06-21-48.mat' ...
    'Results_9306087\res0602_9306087_4_shifts_4_multiQC_50_2025-06-07-07-43.mat'};


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
    % Create dose-volume histogram for target BED
    t = linspace(0, 1.15, n_s);
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
figure
hold on
t = linspace(0, 1.15, n_s);
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
axis([80 135 0 100])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
%legend('JDPO','MIP-JDPO','Location','northeast','FontSize',14);
grid on
%title(strcat('BED DVH in tumor for case ID: ',ptid));
%saveas(gcf,strcat(out_folder,'DTumor_2908_',ptid,'.fig'));
hold off


figure
hold on
t = linspace(0, 1.2, n_s);

try
    % plot(t*100,BED_O{1}(:,1),'LineWidth',2,'LineStyle','--','Color','red','DisplayName','Lung');
    % plot(t*100,BED_O{1}(:,2),'LineWidth',2,'LineStyle','--','Color','blue','HandleVisibility','off');
    L(1:2) = plot(t*100,100*D_O{1}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
    L(3:4) = plot(t*100,100*D_O{1}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
    L(5:6) = plot(t*100,100*D_O{1}(:,3),'g--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
catch
    
end
try
    % plot(t*100,BED_O{2}(:,1),'LineWidth',2,'LineStyle','-','Color','red','DisplayName','Heart');
    % plot(t*100,BED_O{2}(:,2),'LineWidth',2,'LineStyle','-','Color','blue','HandleVisibility','off');
    L(7:8) = plot(t*100,100*D_O{2}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
    L(9:10) = plot(t*100,100*D_O{2}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
    L(11:12) = plot(t*100,100*D_O{2}(:,3),'g--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
catch
    
end
try
    % plot(t*100,BED_O{3}(:,1),'LineWidth',2,'LineStyle',':','Color','red','DisplayName','Esophagus');
    % plot(t*100,BED_O{3}(:,2),'LineWidth',2,'LineStyle',':','Color','blue','HandleVisibility','off');
    L(13:14) = plot(t*100,100*D_O{3}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
    L(15:16) = plot(t*100,100*D_O{3}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
    L(17:18) = plot(t*100,100*D_O{3}(:,3),'g--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
catch
    
end
% for ff = 1:numel(fnames)
%     if ff == 1
%     plot(t,BED_O(:,ff),'LineWidth',2,'LineStyle',':','Color','red');
%     else
%         plot(t,BED_O(:,ff),'LineWidth',2,'Color','red');
%     end
%     %plot(t,BED_O(:,ff),'LineWidth',2);
% end
%axis([0 110 0 80])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
%legend('Lung','Heart','Esophagus','Location','northeast','FontSize',14);
% legend show
%legend(L([2,6]), 'Mandible','Oral','Location','northeast','FontSize',14)
grid on
%title(strcat('BED DVH in lung for case ID: ',ptid));
%saveas(gcf,strcat(out_folder,'DOAR_2908_',ptid,'.fig'));
hold off



