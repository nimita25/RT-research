function y=AmX_robust(x,var)
Dij=var.Dij;
id_obj=var.id_obj;
N_obj=var.N_obj;
N_dij=var.N_Dij;

y=cell(N_obj,N_dij);
for j=1:N_dij
    y0=Dij{j}*double(x);
    for i=1:N_obj
        y{i,j}=y0(id_obj{i,j});
    end
end