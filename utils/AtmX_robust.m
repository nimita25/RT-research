function x=AtmX_robust(y,var)
Dij=var.Dij;
id_obj=var.id_obj;
N_obj=var.N_obj;
w_obj=var.w_obj;
n_c=var.n_c;
c_obj=var.c_obj;
nY=var.nY;
nX=var.nX;
N_Dij=var.N_Dij;

x=zeros([nX,1]);
for j=1:N_Dij
    y0=zeros([nY,1]);
    for i=1:N_obj
%         [i,j]
        y0(id_obj{i,j})=y0(id_obj{i,j})+w_obj(i,j)/n_c(c_obj(i))*y{i,j};
    end
    x=x+Dij{j}'*y0;
end
x=single(x);
