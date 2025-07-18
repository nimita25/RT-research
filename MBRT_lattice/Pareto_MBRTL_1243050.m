
ptid = '1243050';
id = [0 120 240];
all_methods = [4];
nw=2;

% pvdr = zeros(3,2);
% O1max = zeros(3,2);
% 
% load('output_1243050/res_1243050_3_1_2025-07-05-08-08.mat')
% pvdr(1,1) = pvdr3;
% O1max(1,1) = maxD(4); 
% load('output_1243050/res_1243050_3_1_2025-07-05-09-46.mat')
% pvdr(2,1) = pvdr3;
% O1max(2,1) = maxD(4);
% load('output_1243050/res_1243050_3_1_2025-07-05-10-36.mat')
% pvdr(3,1) = pvdr3;
% O1max(3,1) = maxD(4);
% 
% load('output_1243050/res_1243050_3_4_2025-07-05-08-45.mat')
% pvdr(1,2) = pvdr3;
% O1max(1,2) = maxD(4); 
% load('output_1243050/res_1243050_3_4_2025-07-05-09-26.mat')
% pvdr(2,2) = pvdr3;
% O1max(2,2) = maxD(4);
% load('output_1243050/res_1243050_3_4_2025-07-05-10-46.mat')
% pvdr(3,2) = pvdr3;
% O1max(3,2) = maxD(4);
% 
% %[p,p1]=sort(pvdr(:,1),'ascend');
% p1=[1;2;3];


pvdr = zeros(nw,1);
O1max = zeros(nw,1);
O2max = zeros(nw,1);
O1mean = zeros(nw,1);
O2mean = zeros(nw,1);

for mi = 1:numel(all_methods)
    fname = strcat('output_',ptid,'/resPareto_',ptid,'_',num2str(numel(id)),'_',num2str(all_methods(mi)),'.mat');
    load(fname);
    for wi = 1:numel(out_w)
        pvdr(wi,mi) = out_w{wi}.pvdr3;
        O1max(wi,mi) = out_w{wi}.maxD(4);
        O2max(wi,mi) = out_w{wi}.maxD(7);
        O1mean(wi,mi) = out_w{wi}.MD(4);
        O2mean(wi,mi) = out_w{wi}.MD(7);
    end
end

figure
hold on
for mi = 1:numel(all_methods)
plot(pvdr(:,mi),O2mean(:,mi),'-o')
%plot(pvdr(p1,2),O1max(p1,2),'-o')
end
hold off