%load('Results_2286842_enum\resenum0604_2286842_4_shifts_3_50_2025-06-10-00-34.mat')

num_enum = 81;

objs = zeros(num_enum,1);
times = zeros(num_enum,1);
for i = 1:num_enum
    objs(i) = out_enum{i}.obj_total;
    times(i) = out_enum{i}.time;
end

[objs_sorted,ind] = sort(objs);

load('Results_2286842\res0602_2286842_4_shifts_3_QC_50_2025-06-10-13-35.mat')

k = find(objs_sorted == obj_total);


hold on
plot(1:num_enum,objs_sorted*100,'LineWidth',1.5)
plot(k,obj_total*100,'.-','MarkerSize',20,'Color','red')
xlabel('Rank')
ylabel('Obj fn val')
hold off

disp([sum(times), time_mip*10+time_admm]);


