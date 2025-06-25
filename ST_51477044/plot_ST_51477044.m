
RBE = 1.1;
nfrac_list = [5 10 15 20 25 30 35 40 45];
Tl = [7 14 35];
%Tl = [7 ];
Td = [2 20 35];
p1exists = 0;
p2exists = 0;

clear resp1
% Load problem P1 data
resp1 = cell(numel(Tl),numel(Td));
% rhoO, rhoT = 3, 6
% load('.\output\ST-2024-10-09-18-19.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2024-10-09-19-26.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-09-20-33.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2024-10-09-21-40.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2024-10-09-22-47.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2024-10-09-23-53.mat');
% resp1{2,3} = resT;
% load('.\output\ST-2024-10-10-01-00.mat');
% resp1{3,1} = resT;
% load('.\output\ST-2024-10-10-02-06.mat');
% resp1{3,2} = resT;
% load('.\output\ST-2024-10-10-03-12.mat');
% resp1{3,3} = resT;

% % rhoO, rhoT = 3, 4
% load('./output/ST-2024-10-11-14-08.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2024-10-11-15-16.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-11-16-23.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2024-10-11-17-30.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2024-10-11-18-37.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2024-10-11-19-44.mat');
% resp1{2,3} = resT;
% load('.\output\ST-2024-10-11-20-52.mat');
% resp1{3,1} = resT;
% load('.\output\ST-2024-10-11-21-59.mat');
% resp1{3,2} = resT;
% load('.\output\ST-2024-10-11-23-06.mat');
% resp1{3,3} = resT;

% rhoO, rhoT = 6, 4
load('./output/ST-2024-10-12-00-14.mat');
resp1{1,1} = resT;
load('.\output\ST-2024-10-12-01-20.mat');
resp1{1,2} = resT;
load('.\output\ST-2024-10-12-02-27.mat');
resp1{1,3} = resT;
load('.\output\ST-2024-10-12-03-34.mat');
resp1{2,1} = resT;
load('.\output\ST-2024-10-12-04-41.mat');
resp1{2,2} = resT;
load('.\output\ST-2024-10-12-05-48.mat');
resp1{2,3} = resT;
load('.\output\ST-2024-10-12-06-55.mat');
resp1{3,1} = resT;
load('.\output\ST-2024-10-12-08-03.mat');
resp1{3,2} = resT;
load('.\output\ST-2024-10-12-09-10.mat');
resp1{3,3} = resT;
p1exists = 1;

clear resp2
% Load problem P2 data
% load('.\output\ST-2024-10-10-03-53.mat');
% resp2{1,1} = resT;
% load('.\output\ST-2024-10-10-04-35.mat');
% resp2{1,2} = resT;
% load('.\output\ST-2024-10-10-05-17.mat');
% resp2{1,3} = resT;
p2exists = 0;

% Create data for plots for problem P1
if p1exists == 1
BEDt1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
BEDO1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
BEDO11 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
BEDO111 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDt1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDO1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDO11 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDO111 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));

for j = 1:numel(Tl)
for i = 1:numel(Td)
    for t = 1:numel(nfrac_list)
        nfrac = nfrac_list(t);

        % x0 = resp1{j,i}{t}.best_fluencevector;
        % nfrac = resp1{i}{t}.nfrac;
        % x0 = max(0,x0);
        % d = Cost_matrix{1}*x0;
        d = resp1{j,i}{t}.d;
        BED2 = nfrac*(RBE*d + rho(3)*d.^2);
        BED2(c{1}) = BED2(c{1}) - max(nfrac_list(t)-Tl(j)-1,0)*log(2)/Td(i);

        meanBEDtarget = mean(BED2(c{1}));
        meanOAR1 = mean(BED2(c{3}));
        meanOAR2 = mean(BED2(c{4}));
        meanOAR3 = mean(BED2(c{5}));
        
        %BEDO1(t,(j-1)*numel(Tl)+i) = mean(BED2(c{3}))+mean(BED2(c{4}))+mean(BED2(c{5}))+mean(BED2(c{6}))+mean(BED2(c{7}));
        BEDt1(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.output_meanBEDtarget;%meanBEDtarget;
        BEDO1(t,(j-1)*numel(Tl)+i) = meanOAR1;%resp1{j,i}{t}.resd(4);%meanBladder;%resp1{j,i}{t}.resd(5);%sum(resp1{wi}{t}.best_output(4:end));%resp1{wi}{t}.resd(5); 
        BEDO11(t,(j-1)*numel(Tl)+i) = meanOAR2;
        BEDO111(t,(j-1)*numel(Tl)+i) = meanOAR3;
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
        nfrac = nfrac_list(t);

        % x0 = resp1{j,i}{t}.best_fluencevector;
        % nfrac = resp1{i}{t}.nfrac;
        % x0 = max(0,x0);
        % d = Cost_matrix{1}*x0;
        %d = resp2{j,i}{t}.d;
        %BED2 = nfrac*(RBE*d + rho(1)*d.^2);
        %BED2(c{1}) = BED2(c{1}) - max(nfrac_list(t)-Tl(j)-1,0)*log(2)/Td(i);
        
        %meanBEDtarget = mean(BED2(c{1}));
        meanBladder = mean(BED2(c{5}));
        
        BEDt2(t,(j-1)*numel(Tl)+i) = resp2{j,i}{t}.output_meanBEDtarget;%meanBEDtarget;
        BEDO2(t,(j-1)*numel(Tl)+i) = meanBladder;%resp2{j,i}{t}.resd(4);%meanBladder;%resp1{j,i}{t}.resd(5);%sum(resp1{wi}{t}.best_output(4:end));%resp1{wi}{t}.resd(5); 
    end
end
end
end

% Create data for plots for problem P2
% if p2exists == 1
% BEDt2 = zeros(numel(nfrac_list),numel(Td));
% BEDO2 = zeros(numel(nfrac_list),numel(Td));
% 
% for i = 1:numel(Td)
%     for t = 1:numel(nfrac_list)
%         x0 = resp2{i}{t}.best_fluencevector;
%         nfrac = resp2{i}{t}.nfrac;
%         x0 = max(0,x0);
%         d = Cost_matrix{1}*x0;
%         BED2 = nfrac*(RBE*d + rho(1)*d.^2);
%         BED2(c{1}) = BED2(c{1}) - max(nfrac-Tl-1,0)*log(2)/Td(i);
% 
%         meanBEDtarget = mean(BED2(c{1}));
% 
%         BEDt2(t,i) = meanBEDtarget;
%         BEDO2(t,i) = resp2{i}{t}.resd(5);%sum(resp1{wi}{t}.best_output(4:end));%resp1{wi}{t}.resd(5); %5: BED20_bladder
%     end
% end
% end

PlotColors = {'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', "#7E2F8E"};

figure
hold on
legend_name = cell(numel(Tl)*numel(Td),1);
if p1exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
        %plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','Color',PlotColors{i});
        legend_name{(j-1)*numel(Tl)+i} = num2str(Tl(j))+","+num2str(Td(i));
    end
    end
end
if p2exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDt2(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','--','Color',PlotColors{i});
        plot(nfrac_list, BEDt2(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','--');
    end
    end
end
%plot(nfrac_list,resp1{1}{1}.BED_target*ones(numel(nfrac_list),1),'LineStyle',':','Color','black')
xlabel('Number of fractions');
ylabel('mean BED delivered to target');
%legend('Td=2','Td=20','Td=35','Location','southeast')
legend(legend_name,'Location','southeast')
title('Result for prostate case (Tl,Td)')
hold off

figure
hold on
if p1exists == 1
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Marker','.','Color',PlotColors{i});
        plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
    end
    end
    hold off
    figure
    hold on
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Marker','.','Color',PlotColors{i});
        plot(nfrac_list, BEDO11(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
    end
    end
    hold off
    figure
    hold on 
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDO1(:,(j-1)*numel(Tl)+i),'Marker','.','Color',PlotColors{i});
        plot(nfrac_list, BEDO111(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
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
        %plot(nfrac_list, BEDO2(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','--','Color',PlotColors{i});
        plot(nfrac_list, BEDO2(:,(j-1)*numel(Tl)+i),'Marker','.','LineStyle','--');
    end
    end
end
xlabel('Number of fractions');
ylabel('mean BED delivered to femhead');
%legend('Td=2','Td=20','Td=35','Location','southwest')
legend(legend_name,'Location','southwest')
title('Result for prostate case (Tl,Td)')
hold off


% figure
% hold on
% if p1exists == 1
%     for i = 1:numel(Td)
%         plot(nfrac_list, BEDO1(:,i),'Color',PlotColors{i},'Marker','.');
%     end
% end
% if p2exists == 1
%     for i = 1:numel(Td)
%         plot(nfrac_list, BEDO2(:,i),'Color',PlotColors{i},'Marker','.','LineStyle','--');
%     end
% end
% xlabel('Number of fractions');
% ylabel('BED received by OARs');
% legend('Td=2','Td=20','Td=35','Location','northwest')
% title('Result for prostate case with Tl (lag time) = 2')
% hold off


% % Create data for plots for problem P2
% if p2exists == 1

% BEDt2 = zeros(numel(nfrac_list),1);
% BEDO2 = zeros(numel(nfrac_list),1);
% 
% for t = 1:numel(nfrac_list)
%     err = 1e10;
%     for i = 1:numel(resp2)
%         err1 = abs(resp2{i}{t}.BED_target - resp2{i}{t}.output_meanBEDtarget);
%         if err1 < err
%             err = err1;
%             wi = i;
%         end
%     end
%     BEDt2(t) = resp2{wi}{t}.output_meanBEDtarget;
%     BEDO2(t) = resp2{wi}{t}.resd(5); %5: BED20_bladder
% end
% end

% figure
% hold on
% if p1exists == 1
% plot(nfrac_list,BEDt1,'Marker','.');
% plot(nfrac_list,BEDt1b,'Marker','.');
% end
% if p2exists == 1
% plot(nfrac_list,BEDt2,'Marker','.','LineStyle','--');
% end
% plot(nfrac_list,resp1{wi}{t}.BED_target*ones(numel(nfrac_list),1),'LineStyle',':','Color','black')
% xlabel('Number of fractions');
% ylabel('mean BED to target');
% hold off
% 
% figure
% hold on
% if p1exists == 1
% plot(nfrac_list,BEDO1,'Marker','.');
% plot(nfrac_list,BEDO1b,'Marker','.');
% end
% if p2exists == 1
% plot(nfrac_list,BEDO2,'Marker','.','LineStyle','--');
% end
% xlabel('Number of fractions');
% ylabel('BED received by OAR')
% hold off