load('Results_9306087_enum\resenum0604_9306087_4_shifts_3_50_2025-06-06-12-16.mat')

num_enum = 81;

objs = zeros(num_enum,1);
times = zeros(num_enum,1);
for i = 1:num_enum
    objs(i) = out_enum{i}.obj_total;
    times(i) = out_enum{i}.time;
end

[objs_sorted,ind] = sort(objs);

load('Results_9306087\res0602_9306087_4_shifts_3_QC_50_2025-06-06-16-49.mat')
%load('Results_HN02\res0602_HN02_4_shifts_3_QC_50_2025-06-05-08-28.mat')

k = find(objs_sorted == obj_total);


hold on
plot(1:num_enum,objs_sorted*100,'LineWidth',1.5)
plot(k,obj_total*100,'.-','MarkerSize',20,'Color','red')
xlabel('Rank')
ylabel('Obj fn val')
hold off

disp([sum(times), time_mip*10+time_admm]);


