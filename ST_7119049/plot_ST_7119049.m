nfrac_list = [5 10 15 20 25 30 35 40 45 50];
Tl = [7 14 35];
Td = [2 20 35];
p1exists = 0;
p2exists = 0;

clear resp1
% Load problem P1 data
% rhoT = 6, rhoOAR = 3
resp1 = cell(numel(Tl),numel(Td));
% load('.\output\ST-2024-10-06-15-18.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2024-10-07-16-07.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-07-18-51.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2024-10-07-21-34.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2024-10-08-00-17.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2024-10-08-03-00.mat');
% resp1{2,3} = resT;
% load('.\output\ST-2024-10-08-05-43.mat');
% resp1{3,1} = resT;
% load('.\output\ST-2024-10-08-08-27.mat');
% resp1{3,2} = resT;
% load('.\output\ST-2024-10-08-11-11.mat');
% resp1{3,3} = resT;

% Load problem P1 data
% rhoT = 3, rhoOAR = 6
load('.\output\ST-2025-03-17-18-58.mat');
resp1{1,1} = resT;
load('.\output\ST-2025-03-17-21-58.mat');
resp1{1,2} = resT;
load('.\output\ST-2025-03-18-12-37.mat');
resp1{1,3} = resT;
load('.\output\ST-2025-03-18-03-57.mat');
resp1{2,1} = resT;
load('.\output\ST-2025-03-18-06-53.mat');
resp1{2,2} = resT;
load('.\output\ST-2025-03-17-20-27.mat');
resp1{2,3} = resT;
load('.\output\ST-2025-03-17-23-46.mat');
resp1{3,1} = resT;
load('.\output\ST-2025-03-18-03-06.mat');
resp1{3,2} = resT;
load('.\output\ST-2025-03-18-06-25.mat');
resp1{3,3} = resT;
p1exists = 1;

clear resp2
% Load problem P2 data
resp2 = cell(numel(Tl),numel(Td));
% load('.\output\ST-2024-10-06-16-30.mat');
% resp2{1,1} = resT;
% load('.\output\ST-2024-10-06-18-55.mat');
% resp2{1,2} = resT;
% load('.\output\ST-2024-10-06-21-20.mat');
% resp2{1,3} = resT;
% load('.\output\ST-2024-10-06-23-46.mat');
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
BEDO111 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));

for j = 1:numel(Tl)
for i = 1:numel(Td)
    for t = 1:numel(nfrac_list)
        
        BEDt1(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.output_meanBEDtarget;
        BEDO1(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.resd(4);%+resp1{j,i}{t}.resd(6);%resp1{j,i}{t}.resd(4);
        BEDO11(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.resd(6);
        BEDO111(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.resd(8);
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
        
        BEDt2(t,(j-1)*numel(Tl)+i) = resp2{j,i}{t}.output_meanBEDtarget;
        BEDO2(t,(j-1)*numel(Tl)+i) = resp2{j,i}{t}.resd(4);%sum(resp1{wi}{t}.best_output(4:end));%resp1{wi}{t}.resd(5); 
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
        plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
        legend_name{(j-1)*numel(Tl)+i} = num2str(Tl(j))+","+num2str(Td(i));
    end
    end
end

if p2exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        plot(nfrac_list, BEDt2(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
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
%title('Result for lung case (Tl,Td)')
hold off

figure
hold on
if p1exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
        plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
    end
    end
    hold off
    figure
    hold on 
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
        plot(nfrac_list, BEDO11(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
    end
    end
    hold off
end
% if p1exists == 1
%     for j = 1:numel(Tl)
%     for i = 1:numel(Td)
%         %plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
%         plot(nfrac_list, BEDO11(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','-');
%         plot(nfrac_list, BEDO111(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle',':');
%         legend_name{(j-1)*numel(Tl)+i} = num2str(Tl(j))+","+num2str(Td(i));
%     end
%     end
% end

if p2exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        plot(nfrac_list, BEDO2(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
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
%title('Result for lung case (Tl,Td)')
hold off
