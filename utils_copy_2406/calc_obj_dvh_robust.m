function [obj,d]=calc_obj_dvh_robust(x,var)
Dij=var.Dij;
N_obj=var.N_obj;
id_obj=var.id_obj;
w_obj=var.w_obj;
c_obj=var.c_obj;
s_obj=var.s_obj;
N_c=var.N_c;
n_c=var.n_c;
c=var.c;
N_Dij=var.N_Dij;

obj=zeros(N_obj,N_Dij);
d=cell(N_c,N_Dij);
for j=1:N_Dij
    y0=Dij{j}*double(x);
    
    for i=1:N_obj
        obj(i,j)=w_obj(i,j)/n_c(c_obj(i))*sum((y0(id_obj{i,j})-s_obj(i)).^2);
    end
    % obj_total=sum(obj);
    
    for i=1:N_c
        d{i,j}=y0(c{i});
    end
end