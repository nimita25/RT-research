function x=BtX_v3(y,var) % only x-direction total variation
w_i=var.r_i;
w_x=var.r_x;
w_z=var.r_z;
Dij=var.Dij;
intp=var.intp;
tvx=var.tvx;
tvz=var.tvz;
nX=var.nX;
N_dij=var.N_dij;
x=zeros([nX,1]);

for j=1:N_dij
    y1=w_i*y{1,j}/size(intp,1);
    y1=y1+w_x*tvx'*y{2,j}/size(tvx,1);
    y1=y1+w_z*tvz'*y{3,j}/size(tvz,1);
    y0=intp'*y1;
    x=x+Dij{j}'*y0;
end
x=single(x);
end
