function [y,Id]=shrink_gs(x,nnz_x,var)
% vs=struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(N));
%wp=struct('Nx',uint32(nX),'na',uint32(ns),'nx',uint32(n_gs),'isC',uint32(0),'w',ones([nX 1],'single'));
% w=var.w;
ns = var.na;%ns=numel(n_gs);
nX = var.Nx;%nX=size(Dij{1},2);
N = single(var.nx);%n_gs
% N_id = single(var.nx_id);%n_gs_id
% x=single(x);
xs = zeros(ns,1,'single');
m = 0;
for i = 1:ns
    id = m + (1:N(i));
    xs(i) = sum(x(id));% the sum of weights for each energy layer
    m = m + N(i);
end
[tmp,Id] = sort(xs,'descend');
Id = Id(1:nnz_x);
Id_avoid = Id(nnz_x+1:end);
y = zeros(nX,1,'single');
m = 0;
count = 0;
for i = 1:ns
    id = m + (1:N(i));
    if ismember(i,Id)
        y(id) = x(id);
        count = count + 1;
    end
    m = m + N(i);
end
norm(x - y)
end


