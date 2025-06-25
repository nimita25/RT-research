function y=BX_v3(x,var) % total variation on both x&z direction
Dij=var.Dij;
intp=var.intp;
tvx=var.tvx;
tvz=var.tvz;
N_dij=var.N_dij;

y=cell(3,N_dij);
for j=1:N_dij
    y0=Dij{j}*double(x);
    y1=intp*y0;
    y{1,j}=y1;
    y{2,j}=tvx*y1;
    y{3,j}=tvz*y1;
end