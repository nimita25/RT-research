load('Results_7119049_enum/resenum0604_7119049_3_shifts_3_50_2025-06-26-19-27.mat')

num_enum = 27;

objs = zeros(num_enum,1);
times = zeros(num_enum,1);
for i = 1:num_enum
    objs(i) = out_enum{i}.obj_total;
    times(i) = out_enum{i}.time;
end

[objs_sorted,ind] = sort(objs);

load('Results_7119049/res0602_7119049_3_shifts_3_QC_50_2025-06-26-20-30.mat')

k = find(objs_sorted == obj_total);


hold on
plot(1:num_enum,objs_sorted*100,'LineWidth',1.5)
plot(k,obj_total*100,'.-','MarkerSize',20,'Color','red')
xlabel('Rank')
ylabel('Obj fn val')
hold off

disp([sum(times), time_mip*10+time_admm]);


