
ptid='9306087'; 
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=1;
gridDim = '111';
dr = 10;
ddr = 3;
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);
load('9306087_c.mat')

%load('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\output_9306087\res_9306087_4_1_2024-11-13-02-04.mat')

% ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
% body = cst{1,4}{1};
% N_oar = 5;
% oar = cell(N_oar,1);
% oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy; Dmax<54Gy
% oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy; Dmax<36Gy
% oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy; Dmax<54Gy
% oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy; Dmax<54Gy
% oar(5) = cst{44,4}(1); % Brain V12<5cc; Dmax<60Gy

R = [ 4 8 8 ];
fname = strcat('C:\Users\nshinde\Desktop\pMBRT\MBRT_lattice\LVertices\LVertices' ,gridDim, '_',num2str(dr),'_',num2str(ddr),'_9306087_',num2str(R(1)),num2str(R(2)),num2str(R(3)),'.mat');
load(fname);


Vpeak = setdiff(Vpeak,Vlattice{3}{1});
Vvalley = [Vvalley;Vlattice{3}{1}];
c{1} = Vpeak;
c{2} = Vvalley;

px=2;nfrac=10;px0=px*nfrac;pvdr = 5;

%px = px*pvdr;

nm = [3;7];
out_folder = "output/";
resd = cell(3,1);
fnames = {'output_9306087\res_9306087_4_1_2024-11-16-22-37.mat' ...
    'output_9306087\res_9306087_4_1_2024-11-16-22-22.mat',...
    'output_9306087\res_9306087_4_1_2024-11-16-22-52.mat'};%'./output/BED-ADMM-2024-09-02-08-22.mat'};


% Initialize values for DVH plots
n_s = 100;


for ff = 1:numel(fnames)
    
    % Load output here
    fname = fnames{ff};
    disp('Loading output...')
    load(fname);
    
    
    

    % Generating DVH plot for target
    %DVH_eval = BED2(c{1})/BED_target;
    DVH_eval = d(c{1})/(px*pvdr);
    N = numel(DVH_eval);
    % Create dose-volume histogram for target BED
    t = linspace(0, 2, n_s);
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
t = linspace(0, 2, n_s);
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
axis([80 200 0 100])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('pMBRTL-1','pMBRTL-1e','pMBRTL-2','Location','northeast','FontSize',14);
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
    L(11:12) = plot(t*100,100*D_O{2}(:,2),'g--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
catch
    
end

xlabel('BED (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
%legend('Lung','Heart','Esophagus','Location','northeast','FontSize',14);
% legend show
legend(L([2,8]), 'Brainstem','Brain','Location','northeast','FontSize',14)
grid on
%title(strcat('BED DVH in lung for case ID: ',ptid));
%saveas(gcf,strcat(out_folder,'DOAR_2908_',ptid,'.fig'));
hold off



