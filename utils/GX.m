function y=GX(x,var)
% var:wp struct
% y = Wx

% vs=struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(N));
w=var.w;%ones([nX 1]
ns=var.na;%ns=numel(n_gs);
% nX=var.Nx;
N=single(var.nx);%uint32(n_gs)
y=zeros(ns,1,'single');
% x=single(x);

m=0;
for i=1:ns
    id=m+(1:N(i));
    %y(i)=w(id)'*x(id)/N(i);
    y(i)=w(id)'*x(id);
    m=m+N(i);
end
