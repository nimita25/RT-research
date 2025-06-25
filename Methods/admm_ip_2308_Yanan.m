function [xp,s,sg] = admm_ip_2308(ip,var_CG,NNZ,x_init)
%% Parameters set up
mu_min = ip.mu_min;
eps = mu_min/2;
N_dij = numel(var_CG.var_AtA.Dij);
N_obj=var_CG.var_AtA.N_obj;
n_ctv95 = round(var_CG.var_AtA.n_c(1)*0.95);
nX = var_CG.var_AtA.nX;
px = var_CG.var_AtA.s_obj(1);
mup = ip.mup;
var_CG.mu_xs = ip.mup;
update_ac = ip.update_ac;
calc_obj_dvh = ip.calc_obj_dvh;
AtmX = var_CG.AtmX;
AmX = var_CG.AmX;
mu = ip.mu;
N = NNZ;
var_CG.mu_xr = ip.mu;
%Wtx = var_CG.Wtxr;
%Wx = var_CG.Wxr;
wp = var_CG.wpr;
n_gs = wp.nx;
n_el = single(wp.na);
x = x_init;
up = x; zp = up;
%Initialize s
%s = Wx(x,wp);
% [tmp,id] = sort(s,'descend');
% s = zeros(ip.ns,1,'single'); 
% s(id(1:N)) = 1;
s = ones(n_el,1,'single')*N/n_el;
% z = zeros(n_el,1,'single');
% z(1:N) = 1;
z = zeros(n_el,1,'single');
if nargin == 5
    z(isg) = 1;
else
    z(1:N) = 1;
end
u = z; 
sg = cell(N,1);
N_iter = ip.N_iter;

%Save original Dij and use Dij.xj and Dij.sj, where . is Hadamard product
Dij_orig = var_CG.var_AtA.Dij;
size_Dij = size(Dij_orig{1});

%% Optimization 
for n_iter = 1:N_iter
    %% Update Dij according to active energy layers
    nn = 0;
    for ii = 1:n_el
        var_CG.var_AtA.Dij{1}(:,nn+(1:n_gs(ii))) = Dij_orig{1}(:,nn+(1:n_gs(ii)))*double(s(ii));
        nn = nn+n_gs(ii);
    end
    [nY, nX] = size(var_CG.var_AtA.Dij{1});
    var_CG.var_AtA.nX = nX;
    var_CG.var_AtA.nY = nY;

    %% Update active index set 
    [rhs,var_CG.var_AtA] = update_ac(x,var_CG.var_AtA);
    b = AtmX(rhs,var_CG.var_AtA);

    %% Use CG to solve the first subproblem
    [x,cg_err,cg_n] = CG033114(x,b+mup*(zp-up),var_CG.JtJ,var_CG,var_CG.cg_tol,var_CG.cg_iter);
    disp(norm(x));

    %% Update the second subproblem
    xmin = mu_min;
    xt = up + x;
    id = find(xt >= xmin/2);
    zp = zeros(nX,1,'single');
    zp(id) = max(xt(id),xmin);

    %% Calculate D(times)x
    Dx = [];
    nn = 0;
    for ii = 1:n_el
        Dx = [Dx Dij_orig{1}(:,nn+(1:n_gs(ii)))*double(x(nn+(1:n_gs(ii))))];
        nn = nn+n_gs(ii);
    end
    Dx = {Dx};
    var_CG.var_AtA.Dij = Dx;
    [nY, nX] = size(Dx{1});
    var_CG.var_AtA.nX = nX;
    var_CG.var_AtA.nY = nY;

    %% Calculate quantities needed to update s
    [rhs,var_CG.var_AtA]=update_ac(s,var_CG.var_AtA); %%
    bs=AtmX(rhs,var_CG.var_AtA);
    DTD = [];
    for ii = 1:n_el
        tmp_rhs = cell(N_obj,N_dij);
        for j = 1:N_dij
            for i = 1:N_obj
                tmp_rhs{i,j} = Dx{1}(var_CG.var_AtA.id_obj{i,j},ii);
            end
        end
        DTD=[DTD AtmX(tmp_rhs,var_CG.var_AtA)];
    end
    
    %% Update s -> the variable that activates energy layers
    A_lhs = DTD + mu*eye(n_el,n_el);
    b_rhs = bs + mu*(z-u);
    s_temp = CG(double(A_lhs),double(b_rhs),double(s),1e-3,5);
    s = single(s_temp);
    % s = mldivide(A_lhs,b_rhs);

    %% Update auxillary (and binary) variable z
    st = u + s;%
    id = find(st>0.5);
    %[tmp,id] = sort(st,'descend');
    z = zeros(n_el,1,'single'); 
    if ~isempty(id)
        z(id(1:min(N,length(id)))) = 1;
    end

    %% Update the Lagrangian multipliers
    up = up + x - zp;
    u = u + s - z;

    %% Project to MMU constraint
    xp = x;
    xp(xp < eps) = 0;
    xp(intersect(find(xp >= eps),find(xp < mu_min))) = mu_min;

    %% Update Dij according to active energy layers
    nn = 0;
    for ii = 1:n_el
        var_CG.var_AtA.Dij{1}(:,nn+(1:n_gs(ii))) = Dij_orig{1}(:,nn+(1:n_gs(ii)))*double(z(ii));
        nn = nn+n_gs(ii);
    end
    [nY, nX] = size(var_CG.var_AtA.Dij{1});
    var_CG.var_AtA.nX = nX;
    var_CG.var_AtA.nY = nY;

    %% Calculate plan parameter
    %var_CG.var_AtA.Dij = Dij_orig;
    [obj,D] = calc_obj_dvh(xp,var_CG.var_AtA);
    obj_total = sum(obj);  
    D95 = zeros(N_dij,1);
    Dmax = zeros(N_dij,1);
    for i = 1:N_dij
        y2 = sort(D{1,i},'descend');
        D95(i) = y2(n_ctv95)/px;
        Dmax(i) = y2(1)/px;
    end
   %% Display
   disp(num2str([n_iter mean(obj_total) mean(D95) mean(Dmax)]))

   %% Find beam intensities only for active energy layers (others=0)
   if n_iter == N_iter
       sg = find(z==1);
       nn = 0;
       for ii = 1:n_el
           xp(nn+(1:n_gs(ii))) = xp(nn+(1:n_gs(ii)))*z(ii);
           nn = nn+n_gs(ii);
       end
       %[xp,Id] = shrink_gs(xp,N,wp);
    end

end





