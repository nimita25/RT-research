clc;clear;close all;
ptid='2286842'; % (!)
id=[45 135 225 315];
px=2.12;nfrac=33;

name1=['res_' ptid '_4_1']; % convention
if 1
name2=['res_' ptid '_4_5']; % convention
name3=['res_' ptid '_4_4_6_6_0_0.1']; % convention
name4=['res_' ptid '_4_7']; % convention
else
name2=['res_' ptid '_4_7']; % convention
name3=['res_' ptid '_4_6_4_4_0_0.1']; % convention
name4=['res_' ptid '_4_5']; % convention
end

load(['D:\KUMC\newdata101921\' ptid '\' ptid '.mat'],'ct','cst');
load([ptid '_c.mat'],'c');
load([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat']);
load([ptid '\' ptid '_' num2str(numel(id)) '_tvm.mat']);
nc=numel(c);

res=zeros(9,20);
%%%%%%%%%% plan 1 %%%%%%%%%%
load([ptid '\' name1 '.mat']);

d=d/px;
d_oar1=cell(nc,1);
for j=1:nc
    d_oar1{j}=d(c{j});
end
di11=reshape(intp_45*d(:),slice_Dim);
di12=reshape(intp_135*d(:),slice_Dim);
di13=reshape(intp_225*d(:),slice_Dim);
di14=reshape(intp_315*d(:),slice_Dim);
factor1=factor;
[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,Dpeak1,Dvalley1,PVDR1,SF1,Dpeak2,Dvalley2,PVDR2,SF2,...
    Dpeak3,Dvalley3,PVDR3,SF3,Dpeak4,Dvalley4,PVDR4,SF4,Delivery_Time]=calcpara_2286842_4(d,px,nfrac,c,di11,di12,di13,di14,x);
res(1,:)=[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,PVDR1,PVDR2,...
    PVDR3,PVDR4,SF1,SF2,SF3,SF4,Delivery_Time];

%%%%%%%%%% plan 2 %%%%%%%%%%
load([ptid '\' name2 '.mat']);
figure;imshow3D(d,[0,1]);
d=d/px;
d_oar2=cell(nc,1);
for j=1:nc
    d_oar2{j}=d(c{j});
end
di21=reshape(intp_45*d(:),slice_Dim);
di22=reshape(intp_135*d(:),slice_Dim);
di23=reshape(intp_225*d(:),slice_Dim);
di24=reshape(intp_315*d(:),slice_Dim);
factor2=factor;
[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,Dpeak1,Dvalley1,PVDR1,SF1,Dpeak2,Dvalley2,PVDR2,SF2,...
    Dpeak3,Dvalley3,PVDR3,SF3,Dpeak4,Dvalley4,PVDR4,SF4,Delivery_Time]=calcpara_2286842_4(d,px,nfrac,c,di21,di22,di23,di24,x);
res(2,:)=[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,PVDR1,PVDR2,...
    PVDR3,PVDR4,SF1,SF2,SF3,SF4,Delivery_Time];

%%%%%%%%%% plan 3 %%%%%%%%%%
load([ptid '\' name3 '.mat']);
% figure;imshow3D(d,[0,1]);
d=d/px;
d_oar3=cell(nc,1);
for j=1:nc
    d_oar3{j}=d(c{j});
end
di31=reshape(intp_45*d(:),slice_Dim);
di32=reshape(intp_135*d(:),slice_Dim);
di33=reshape(intp_225*d(:),slice_Dim);
di34=reshape(intp_315*d(:),slice_Dim);
factor3=factor;
[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,Dpeak1,Dvalley1,PVDR1,SF1,Dpeak2,Dvalley2,PVDR2,SF2,...
    Dpeak3,Dvalley3,PVDR3,SF3,Dpeak4,Dvalley4,PVDR4,SF4,Delivery_Time]=calcpara_2286842_4(d,px,nfrac,c,di31,di32,di33,di34,x);
res(3,:)=[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,PVDR1,PVDR2,...
    PVDR3,PVDR4,SF1,SF2,SF3,SF4,Delivery_Time];

%%%%%%%%%% plan 4 %%%%%%%%%%
load([ptid '\' name4 '.mat']);
% figure;imshow3D(d,[0,1]);
d=d/px;
d_oar4=cell(nc,1);
for j=1:nc
    d_oar4{j}=d(c{j});
end
di41=reshape(intp_45*d(:),slice_Dim);
di42=reshape(intp_135*d(:),slice_Dim);
di43=reshape(intp_225*d(:),slice_Dim);
di44=reshape(intp_315*d(:),slice_Dim);
factor4=factor;
[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,Dpeak1,Dvalley1,PVDR1,SF1,Dpeak2,Dvalley2,PVDR2,SF2,...
    Dpeak3,Dvalley3,PVDR3,SF3,Dpeak4,Dvalley4,PVDR4,SF4,Delivery_Time]=calcpara_2286842_4(d,px,nfrac,c,di41,di42,di43,di44,x);
res(4,:)=[D98,Dmax,CI,Dmean_body,Dmean_oar1,Dmean_oar2,Dmean_oar3,Dmean_oar4,Dmean_oar5,Dmean_oar6,Dmean_oar7,PVDR1,PVDR2,...
    PVDR3,PVDR4,SF1,SF2,SF3,SF4,Delivery_Time];

idx=1:slice_Dim(1);
idz=75:95;
if 1
figure;imshow(cat(1,cat(2,di11(:,idz),di21(:,idz),di31(:,idz)),...
    cat(2,di12(:,idz),di22(:,idz),di32(:,idz)),...
    cat(2,di12(:,idz),di22(:,idz),di32(:,idz)),...
    cat(2,di14(:,idz),di24(:,idz),di34(:,idz))),[0 0.6]);
% figure;imshow(cat(2,di11(:,idz),di21(:,idz),di31(:,idz)),[0 0.8]);
% figure;imshow(cat(2,di12(:,idz),di22(:,idz),di32(:,idz)),[0 0.6]);
% figure;imshow(cat(2,di13(:,idz),di23(:,idz),di33(:,idz)),[0 0.4]);
end


if 1
idk=85;
figure;
hold on;
plot(idx,di31(:,idk),'r',idx,di21(:,idk),'b',idx,di11(:,idk),'k','linewidth',2)
hold off;
xlim([0 240])
ylim([0 0.7])
yticks([0 0.2 0.4 0.6 0.8])
yticklabels({'0','20','40','60','80'})
xlabel('x-position') 
ylabel('Dose (%)') 
set(gca,'FontSize', 48);
grid on
box on;
grid minor
set(0,'defaultfigurecolor','w')

figure;
hold on;
plot(idx,di32(:,idk),'r',idx,di22(:,idk),'b',idx,di12(:,idk),'k','linewidth',2)
hold off;
xlim([0 240])
ylim([0 0.7])
yticks([0 0.2 0.4 0.6 0.8])
yticklabels({'0','20','40','60','80'})
xlabel('x-position') 
ylabel('Dose (%)') 
set(gca,'FontSize', 48);
grid on
box on;
grid minor
set(0,'defaultfigurecolor','w')

figure;
hold on;
plot(idx,di13(:,idk),'k',idx,di23(:,idk),'b',idx,di33(:,idk),'r','linewidth',2)
hold off;
xlim([0 240])
ylim([0 0.7])
yticks([0 0.2 0.4 0.6 0.8])
yticklabels({'0','20','40','60','80'})
xlabel('x-position') 
ylabel('Dose (%)') 
set(gca,'FontSize', 48);
grid on
box on;
grid minor
set(0,'defaultfigurecolor','w')

figure;
hold on;
plot(idx,di34(:,idk),'r',idx,di24(:,idk),'b',idx,di14(:,idk),'k','linewidth',2)
hold off;
xlim([0 240])
ylim([0 0.7])
yticks([0 0.2 0.4 0.6 0.8])
yticklabels({'0','20','40','60','80'})
xlabel('x-position') 
ylabel('Dose (%)') 
set(gca,'FontSize', 48);
grid on
box on;
grid minor
set(0,'defaultfigurecolor','w')
end

if 1
n=100;
t=(0:1/(n-1):1)*1.2;

for j=[7,8]
figure;hold on;
set(0,'defaultfigurecolor','w') 

    oar1=zeros(n,1);
    oar2=zeros(n,1);
    oar3=zeros(n,1);
    oar4=zeros(n,1);
    oar5=zeros(n,1);
    n_oar=numel(c{j});
    for i=1:n
        oar1(i)=numel(find(d_oar1{j}>=t(i)))/n_oar;
        oar2(i)=numel(find(d_oar2{j}>=t(i)))/n_oar;
        oar3(i)=numel(find(d_oar3{j}>=t(i)))/n_oar;
%         oar4(i)=numel(find(d_oar4{j}>=t(i)))/n_oar;
%         oar5(i)=numel(find(d_oar5{j}>=t(i)))/n_oar;
    end
    
    plot(t,oar1,'k','LineWidth',3);
    plot(t,oar2,'b','LineWidth',3);
    plot(t,oar3,'r','LineWidth',3);
%     plot(t,oar4,'y','LineWidth',3);

hold off
set(gca,'FontSize', 24);
grid on
grid minor

if j==1
xlim([0 1.2]);
set(gca,'XTick',0:0.3:1.2)
set(gca,'XTickLabel',0:30:120);
ylim([0 1]);
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabel',0:20:100);
elseif j==8
xlim([0 0.5]);
set(gca,'XTick',0:0.1:1)
set(gca,'XTickLabel',0:10:100);
ylim([0 0.4]);
set(gca,'YTick',0:0.1:0.5)
set(gca,'YTickLabel',0:10:50);
elseif j==7
xlim([0 0.7]);
set(gca,'XTick',0:0.1:1)
set(gca,'XTickLabel',0:10:100);
ylim([0 0.5]);
set(gca,'YTick',0:0.1:0.5)
set(gca,'YTickLabel',0:10:50);
end

set(0,'defaultfigurecolor','w') 
xlabel('Dose (%)');
ylabel('Volume (%)');
hold off;

end
end