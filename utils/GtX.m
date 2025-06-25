function x=GtX(y,var)
% vs=struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(N));
% y=single(y);
w=var.w;
ns=var.na;
nX=var.Nx;
N=single(var.nx);
x=zeros(nX,1,'single');


m=0;
for i=1:ns
    id=m+(1:N(i));
    %x(id)=w(id)*y(i)/N(i);
    x(id)=w(id)*y(i);
    m=m+N(i);
end
