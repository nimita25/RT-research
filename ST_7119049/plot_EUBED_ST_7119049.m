nfrac_list = [5 10 15 20 25 30 35 40 45 50];
Tl = [7 14 35];
Td = [2 20 35];
p1exists = 0;
p2exists = 0;

ptid = '7119049';
folder = ['../' ptid '/'];
%f = functionsContainer;
load([folder ptid '.mat'], 'ct', 'cst');
ctv1 = cst{9,4}{1};   % ptv60, 2*30
body = cst{1,4}{1};
N_oar = 3;
oar = cell(N_oar,1);
oar(1)=cst{7,4}(1); % lung Dmean<18Gy, V20<30% (D30<12Gy)
oar(2)=cst{4,4}(1); % heart V45<60% (D60<27Gy)
oar(3)=cst{3,4}(1); % esophagus Dmean<20
n_rho = numel(oar)+2;
rho = (1/3)*ones([n_rho,1]);
rho(1) = 1/6;
alpha = 0.04;
rho_target = rho(1);
RBE = 1.1;
% Calculate BED for target (rhs of equality constraint)
t_px = 60; %total prescription dose
nfrac = 30; % number of fraction
px = t_px/nfrac; % prescription dose
ctv = cell(1,1);
ctv{1} = ctv1; %row indices of Dij corresponding to target/tumor
n_oar = zeros(N_oar, 1);
for i = 1:N_oar
    oar{i} = setdiff(oar{i}, ctv1);  %row indices of Dij corresponding to OAR
    n_oar(i) = numel(oar{i}); %number of voxels in each OAR
end
c = [ctv; {body}; oar;]; %cell containing row indices corresponding to OAR, tumor, body

n_c = zeros([numel(c) 1]); %number of voxels in each DVH-max, DVH-mean OAR
for i = 1:numel(n_c)
    n_c(i) = numel(c{i});
end

clear resp1
% Load problem P1 data
% rhoT = 6, rhoOAR = 3
resp1 = cell(numel(Tl),numel(Td));
load('.\output\ST-2024-10-06-15-18.mat');
resp1{1,1} = resT;
load('.\output\ST-2024-10-07-16-07.mat');
resp1{1,2} = resT;
load('.\output\ST-2024-10-07-18-51.mat');
resp1{1,3} = resT;
load('.\output\ST-2024-10-07-21-34.mat');
resp1{2,1} = resT;
load('.\output\ST-2024-10-08-00-17.mat');
resp1{2,2} = resT;
load('.\output\ST-2024-10-08-03-00.mat');
resp1{2,3} = resT;
load('.\output\ST-2024-10-08-05-43.mat');
resp1{3,1} = resT;
load('.\output\ST-2024-10-08-08-27.mat');
resp1{3,2} = resT;
load('.\output\ST-2024-10-08-11-11.mat');
resp1{3,3} = resT;

% Load problem P1 data
% rhoT = 3, rhoOAR = 6
% load('.\output\ST-2025-03-17-18-58.mat');
% resp1{1,1} = resT;
% load('.\output\ST-2025-03-17-21-58.mat');
% resp1{1,2} = resT;
% load('.\output\ST-2024-10-07-18-51.mat');
% resp1{1,3} = resT;
% load('.\output\ST-2025-03-18-03-57.mat');
% resp1{2,1} = resT;
% load('.\output\ST-2025-03-18-06-53.mat');
% resp1{2,2} = resT;
% load('.\output\ST-2025-03-17-20-27.mat');
% resp1{2,3} = resT;
% load('.\output\ST-2025-03-17-23-46.mat');
% resp1{3,1} = resT;
% load('.\output\ST-2025-03-18-03-06.mat');
% resp1{3,2} = resT;
% load('.\output\ST-2025-03-18-06-25.mat');
% resp1{3,3} = resT;
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
EUBEDt1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDO1 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDO11 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));
EUBEDO111 = zeros(numel(nfrac_list),numel(Td)*numel(Tl));

for j = 1:numel(Tl)
for i = 1:numel(Td)
    for t = 1:numel(nfrac_list)
        nfrac = nfrac_list(t);
        d = resp1{j,i}{t}.d;
        BED2 = nfrac*(RBE*d + rho(1)*d.^2);
        BED2(c{1}) = BED2(c{1}) - max(nfrac_list(t)-Tl(j)-1,0)*log(2)/Td(i);
        
        BEDt1(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.output_meanBEDtarget;
        BEDO1(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.resd(4);%+resp1{j,i}{t}.resd(6);%resp1{j,i}{t}.resd(4);
        BEDO11(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.resd(6);
        BEDO111(t,(j-1)*numel(Tl)+i) = resp1{j,i}{t}.resd(8);

        for nc = [1 ]
        n_s = 100;
        DVHeval = BED2(c{nc})/max(BED2(c{nc}));
        N = numel(DVHeval);
        tt = linspace(0,1,n_s);
        v = zeros(n_s,1);
        BEDtmp = zeros(n_s,1);
        for ii = 1:n_s
            if ii<n_s
            v(ii) = numel(find((DVHeval>=tt(ii)& DVHeval<=tt(ii+1))))/N;
            else
            v(ii) = numel(find((DVHeval>=tt(ii))))/N;
            end
            BEDtmp(ii) = tt(ii)*max(BED2(c{nc}));
            if nc == 1
                EUBEDt1(t,(j-1)*numel(Tl)+i) = -(1/alpha)*log(sum(exp(-alpha*BEDtmp).*v));
            elseif nc == 3
                EUBEDO1(t,(j-1)*numel(Tl)+i) = -(1/alpha)*log(sum(exp(-alpha*BEDtmp).*v));
            elseif nc == 4
                EUBEDO11(t,(j-1)*numel(Tl)+i) = -(1/alpha)*log(sum(exp(-alpha*BEDtmp).*v));
            else
                EUBEDO111(t,(j-1)*numel(Tl)+i) = -(1/alpha)*log(sum(exp(-alpha*BEDtmp).*v));
            end
            
        end
        end
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
        plot(nfrac_list, EUBEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
        %plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2,'LineStyle',':');
        legend_name{(j-1)*numel(Tl)+i} = num2str(Tl(j))+","+num2str(Td(i));
    end
    end
    for j = 1:numel(Tl)
    for i = 1:numel(Td)
        %plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Color',PlotColors{i},'Marker','.');
        %plot(nfrac_list, EUBEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2);
        plot(nfrac_list, BEDt1(:,(j-1)*numel(Tl)+i),'Marker','.','LineWidth',2,'LineStyle',':');
        legend_name{(j-1)*numel(Tl)+i} = num2str(Tl(j))+","+num2str(Td(i));
    end
    end
end
plot(nfrac_list,resp1{1}{1}.BED_target*ones(numel(nfrac_list),1),'LineStyle',':','Color','black')


hold off

