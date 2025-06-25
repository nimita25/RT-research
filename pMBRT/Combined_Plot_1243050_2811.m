
ptid='1243050'; 
addpath('C:\Users\nshinde\Desktop\pMBRT\pMBRT');
method=1;
gridDim = '113';
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid' gridDim '.mat']);
load([ptid '_c.mat'],'c');


px=6;%nfrac=10;px0=px*nfrac;pvdr = 5;



nm = [3;5];
out_folder = "output/";
resd = cell(3,1);
fnames = {'../1243050/res_1243050_3_3_2_ctc337.mat' ...
    '../1243050/res_1243050_3_3_2_ctc333.mat'};%'./output/BED-ADMM-2024-09-02-08-22.mat'};


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
    t = linspace(0, 1.3, n_s);
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
t = linspace(0, 1.3, n_s);
tid = 100;
for ff = 1:numel(fnames)
    if ff == 1
    plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','red');
    else
        plot(t(1:tid)*100,100*D_T((1:tid),ff),'LineWidth',2,'Color','blue');
    end
end
axis([80 125 0 100])
xlabel('Dose (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
legend('JDPO','MIP-JDPO','Location','northeast','FontSize',14);
grid on
%title(strcat('BED DVH in tumor for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DTumor_2908_',ptid,'.fig'));
hold off


figure
hold on
t = linspace(0, 1.2, n_s);

try
    % plot(t*100,BED_O{1}(:,1),'LineWidth',2,'LineStyle','--','Color','red','DisplayName','Lung');
    % plot(t*100,BED_O{1}(:,2),'LineWidth',2,'LineStyle','--','Color','blue','HandleVisibility','off');
    L(1:2) = plot(t*100,100*D_O{1}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
    L(3:4) = plot(t*100,100*D_O{1}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle','-');
catch
    
end
try
    % plot(t*100,BED_O{2}(:,1),'LineWidth',2,'LineStyle','-','Color','red','DisplayName','Heart');
    % plot(t*100,BED_O{2}(:,2),'LineWidth',2,'LineStyle','-','Color','blue','HandleVisibility','off');
    L(5:6) = plot(t*100,100*D_O{2}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
    L(7:8) = plot(t*100,100*D_O{2}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle','--');
catch
    
end
try
    % plot(t*100,BED_O{3}(:,1),'LineWidth',2,'LineStyle',':','Color','red','DisplayName','Esophagus');
    % plot(t*100,BED_O{3}(:,2),'LineWidth',2,'LineStyle',':','Color','blue','HandleVisibility','off');
    L(9:10) = plot(t*100,100*D_O{3}(:,1),'r--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
    L(11:12) = plot(t*100,100*D_O{3}(:,2),'b--',nan,nan,'k--','LineWidth',2,'LineStyle',':');
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
xlabel('BED (%)','FontSize',18);
ylabel('Volume (%)','FontSize',18);
%legend('Lung','Heart','Esophagus','Location','northeast','FontSize',14);
% legend show
legend(L([2,6]), 'Large Bowel','Spinal Cord','Location','northeast','FontSize',14)
grid on
%title(strcat('BED DVH in lung for case ID: ',ptid));
saveas(gcf,strcat(out_folder,'DOAR_2908_',ptid,'.fig'));
hold off



