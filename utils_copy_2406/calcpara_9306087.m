function [D95,Dmax,CI,Dmax_oar1,Dmean_oar1,V10_oar1,Dmax_oar2,Dmax_oar3,Dmax_oar4,V12_oar5,Dmean_body]=calcpara_9306087(nfrac,d,px,ctv,oar1,oar2,oar3,oar4,oar5,body)

% nfrac=4;
d=d(:);
px = 2;
d2=d(ctv);
Dmax=max(d2)/px*100;

d3=sort(d2,'descend');
D95=d3(round(0.95*numel(ctv)));

CI=sum(d2>=px)^2/(sum(d>=px)*numel(ctv));

d2=d(oar1);
% Dmax_oar1=max(d2)*nfrac;
Dmax_oar1 = max(d2)/px;
Dmean_oar1 = mean(d2)/px;
V10_oar1=sum(d2>=10/nfrac)*0.3^3;
% Dmean_oar1=mean(d2)*nfrac;
% d3=sort(d2,'descend');
% Dmean_oar1=d3(round(0.1/0.3^3))*nfrac;

d2=d(oar2);
Dmax_oar2=max(d2)*nfrac;

d2=d(oar3);
Dmax_oar3=max(d2)*nfrac;

d2=d(oar4);
Dmax_oar4=max(d2)*nfrac;

d2=d(oar5);
V12_oar5=sum(d2>=12/nfrac)*0.3^3;

d2=d(body);
% Dmean_body=mean(d2)*nfrac;
Dmean_body = mean(d2)/px;

