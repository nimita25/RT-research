nfrac_list = [5 10 15 20 25 30 35 40 45 50 55 60];
Tl = [7 14 35];
Td = [2 20 35];
p1exists = 0;
p2exists = 0;

clear resp1
% Load problem P1 data
resp1 = cell(numel(Tl),numel(Td));
% rhoO, rhoT = 3, 10
load('.\output\ST-2024-10-08-15-30.mat');
resp1{1,1} = resT;
load('.\output\ST-2024-10-08-17-50.mat');
resp1{1,2} = resT;
load('.\output\ST-2024-10-08-20-10.mat');
resp1{1,3} = resT;
load('.\output\ST-2024-10-08-22-30.mat');
resp1{2,1} = resT;
load('.\output\ST-2024-10-09-10-56.mat');
resp1{2,2} = resT;
load('.\output\ST-2024-10-09-12-30.mat');
resp1{2,3} = resT;
load('.\output\ST-2024-10-09-13-57.mat');
resp1{3,1} = resT;
load('.\output\ST-2024-10-09-15-35.mat');
resp1{3,2} = resT;
load('.\output\ST-2024-10-09-17-02.mat');
resp1{3,3} = resT;

% % rhoO, rhoT = 6, 8
% load('.\output\ST-2024-10-12-17-20.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2024-10-12-18-47.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-12-20-14.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2024-10-12-21-41.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2024-10-12-23-08.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2024-10-13-00-35.mat');
% resp1{2,3} = resT;
% load('.\output\ST-2024-10-13-02-02.mat');
% resp1{3,1} = resT;
% load('.\output\ST-2024-10-13-03-29.mat');
% resp1{3,2} = resT;
% load('.\output\ST-2024-10-13-04-56.mat');
% resp1{3,3} = resT;

% rhoO, rhoT = 3, 8
% load('.\output\ST-2024-10-13-06-23.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2024-10-13-07-49.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-13-20-05.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2024-10-13-21-32.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2024-10-13-22-59.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2024-10-14-00-26.mat');
% resp1{2,3} = resT;
% load('.\output\ST-2024-10-14-01-53.mat');
% resp1{3,1} = resT;
% load('.\output\ST-2024-10-14-03-20.mat');
% resp1{3,2} = resT;
% load('.\output\ST-2024-10-14-04-46.mat');
% resp1{3,3} = resT;
% nfrac_list = [30];
% load('.\output\ST-2024-10-14-08-57.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2024-10-14-09-05.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-14-09-13.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2024-10-14-09-20.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2024-10-14-09-28.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2024-10-14-09-36.mat');
% resp1{2,3} = resT;
p1exists = 1;

clear resp2
% Load problem P2 data
resp2 = cell(numel(Tl),numel(Td));
% load('.\output\ST-2024-10-08-16-21.mat');
% resp2{1,1} = resT;
% load('.\output\ST-2024-10-08-18-42.mat');
% resp2{1,2} = resT;
% load('.\output\ST-2024-10-08-21-02.mat');
% resp2{1,3} = resT;
% load('.\output\ST-2024-10-08-23-22.mat');
% resp2{2,1} = resT;
% load('.\output\ST-2024-10-07-02-12.mat');
% resp2{2,2} = resT;
% load('.\output\ST-2024-10-07-04-37.mat');
% resp2{2,3} = resT;
% load('.\output\ST-2024-10-07-07-02.mat');
% resp2{3,1} = resT;
% load('.\output\ST-2024-10-07-09-26.mat');
% resp2{3,2} = resT;
% load('.\output\ST-2024-10-07-11-51.mat');
% resp2{3,3} = resT;
p2exists = 0;

% Create data for plots for problem P1
if p1exists == 1
BEDt1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
BEDO1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
BEDO11 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));

for j = 1:numel(Tl)
for i = 1:numel(Td)
    for t = 1:numel(nfrac_list)
        
        BEDt1(t,(j-1)*numel(Td)+i) = resp1{j,i}{t}.output_meanBEDtarget;
        BEDO1(t,(j-1)*numel(Td)+i) = resp1{j,i}{t}.resd(5);%sum(resp1{wi}{t}.best_output(4:end));%resp1{wi}{t}.resd(5); 
        BEDO11(t,(j-1)*numel(Td)+i) = resp1{j,i}{t}.resd(11);
    end
end
end
end



% Create data for plots for problem P2
if p2exists == 1
BEDt2 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
BEDO2 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));

for j = 1:numel(Tl)
for i = 1:numel(Td)
    for t = 1:numel(nfrac_list)
        
        BEDt2(t,(j-1)*numel(Td)+i) = resp2{j,i}{t}.output_meanBEDtarget;
        BEDO2(t,(j-1)*numel(Td)+i) = resp2{j,i}{t}.resd(5);%sum(resp1{wi}{t}.best_output(4:end));%resp1{wi}{t}.resd(5); 
    end
end
end
end

PlotColors = {'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', "#7E2F8E"};

figure
hold on
legend_name = cell(numel(Tl)*numel(Td),1);
if p1exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
        plot(nfrac_list, BEDt1(:,(j-1)*numel(Td)+i),'Marker','.','LineWidth',2);
        legend_name{(j-1)*numel(Td)+i} = num2str(Tl(j))+","+num2str(Td(i));
    end
    end
end

if p2exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        plot(nfrac_list, BEDt2(:,(j-1)*numel(Td)+i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
        %plot(nfrac_list, BEDt2(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','--');
    end
    end
end
% if p2exists == 1
%     for i = 1:numel(Td)
%         plot(nfrac_list, BEDt2(:,i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
%     end
% end
plot(nfrac_list,resp1{1}{1}.BED_target*ones(numel(nfrac_list),1),'LineStyle',':','Color','black')
%xlabel('Number of fractions');
%ylabel('mean BED delivered to target');
%legend('Td=2','Td=20','Td=35','Location','southeast')
%legend(legend_name,'Location','southeast')
%title('Result for 9306087 case: a/b = 8 (target), 6 (OAR) (Tl,Td)')
hold off

figure
hold on
if p1exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
        plot(nfrac_list, BEDO1(:,(j-1)*numel(Td)+i),'Marker','.','LineWidth',2);
    end
    end
end
hold off
figure
hold on
if p1exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
        plot(nfrac_list, BEDO11(:,(j-1)*numel(Td)+i),'Marker','.','LineStyle','-','LineWidth',2);
    end
    end
end
if p2exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        plot(nfrac_list, BEDO2(:,(j-1)*numel(Td)+i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
        %plot(nfrac_list, BEDO2(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','--');
    end
    end
end
% if p2exists == 1
%     for i = 1:numel(Td)
%         plot(nfrac_list, BEDO2(:,i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
%     end
% end
%xlabel('Number of fractions');
%ylabel('mean BED delivered to OAR');
%legend('Td=2','Td=20','Td=35','Location','southwest')
%legend(legend_name,'Location','southwest')
%title('Result for 9306087 case: a/b = 8 (target), 6 (OAR) (Tl,Td)')
hold off
