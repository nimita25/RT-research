load('Results_HN02_enum\resenum0604_HN02_4_shifts_3_50_2025-06-05-01-10.mat')

num_enum = 81;

objs = zeros(num_enum,1);
times = zeros(num_enum,1);
for i = 1:num_enum
    objs(i) = out_enum{i}.obj_total;
    times(i) = out_enum{i}.time;
end

[objs_sorted,ind] = sort(objs);

load('Results_HN02\res0602_HN02_4_shifts_3_QC_50_2025-06-05-09-08.mat')
%load('Results_HN02\res0602_HN02_4_shifts_3_QC_50_2025-06-05-08-28.mat')

k = find(objs_sorted == obj_total);


hold on
plot(1:num_enum,objs_sorted,'LineWidth',1.5)
plot(k,obj_total,'.-','MarkerSize',20,'Color','red')
xlabel('Rank')
ylabel('Obj fn val')
hold off

disp([sum(times), time_mip*10+time_admm]);


