function [dmf,dmf_dr,dmf_d]=flash_dmf_logistic(dr,d,var)

ctv=var.ctv;
d_c=var.d_c;
dr_c=var.dr_c;
dmf0=var.dmf0;
k=var.k;

% close all;
% k=10;
% c=1.5;
% % c=2;
% d_c=8;
% dr_c=40;

kd=k/d_c;
kdr=k/dr_c;

% n=1000;
% d=d_c*((0:1/n:2)-1);
% % d=d_c;
% dr=dr_c*((0:1/n:2)-1);
% dr=dr_c;

% d=max(d,0);
% dr=max(dr,0);


dmf=1+(1/dmf0-1)*(1+exp(-kd*d)).^(-0.5).*(1+exp(-kdr*dr)).^(-0.5);

dmf(ctv)=1;

dmf_d=(1/dmf0-1)*(0.5*kd*exp(-kd*d).*(1+exp(-kd*d)).^(-1.5)).*(1+exp(-kdr*dr)).^(-0.5);
dmf_dr=(1/dmf0-1)*(1+exp(-kd*d)).^(-0.5).*(0.5*kdr*exp(-kdr*dr).*(1+exp(-kdr*dr)).^(-1.5));

dmf_d(ctv)=0;
dmf_dr(ctv)=0;

dmf(~isfinite(dmf))=1;
dmf_d(~isfinite(dmf_d))=0;
dmf_dr(~isfinite(dmf_dr))=0;


% if sum(isnan(dmf))>0
%     disp(['nan:' num2str(sum(isnan(dmf)))]);
%     dmf(isnan(dmf))=1;
%     dmf_d(isnan(dmf))=0;
%     dmf_dr(isnan(dmf))=0;
% end
% if sum(isinf(dmf))>0
%     disp(['inf:' num2str(sum(isinf(dmf)))]);
%     dmf(isinf(dmf))=1;
%     dmf_d(isinf(dmf))=0;
%     dmf_dr(isinf(dmf))=0;
% end



% [sum(isnan(dmf)) sum(isnan(dmf_d)) sum(isnan(dmf_dr))]

% nY=numel(d);
% dmf=ones(nY,1);
% dmf_d=zeros(nY,1);
% dmf_dr=zeros(nY,1);

% figure;plot(dr,dmf);
% 
% dmf_d2=(dmf(2:end)-dmf(1:end-1))/(d(2)-d(1));
% figure;plot(dmf_d(1:end-1));
% figure;plot(dmf_d2);

% dmf_dr2=(dmf(2:end)-dmf(1:end-1))/(dr(2)-dr(1));
% figure;plot(dmf_dr(1:end-1));
% figure;plot(dmf_dr2);


