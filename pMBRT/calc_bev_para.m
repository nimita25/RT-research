function [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px)

D_bev = mean(dk)*100/px;
dk1=sort(dk,'descend');
dk2=dk1(dk1>0.01);
d10 = ceil(0.1*numel(dk2));
d80 = ceil(0.8*numel(dk2));
pvdr1 = dk2(d10)/dk2(d80);
d20 = ceil(0.2*numel(dk2));
d70 = ceil(0.7*numel(dk2));

pvdr2 = mean(dk2(d10:d20))/mean(dk2(d70:d80));