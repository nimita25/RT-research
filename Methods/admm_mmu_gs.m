function [xp,sg] = admm_mmu_gs(ip,var_CG,NNZ,x_init)
%% Parameters set up
mu_min = ip.mu_min;
eps = mu_min/2;
N_dij = numel(var_CG.var_AtA.Dij);
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
Wtx = var_CG.Wtxr;
Wx = var_CG.Wxr;
wp = var_CG.wpr;
x = x_init;
up = x; zp = up;
u = Wx(x,wp); z = u;
sg = cell(N,1);
n_iter = 0;
N_iter = ip.N_iter;

%% Optimization
disp(size(var_CG.var_AtA.Dij));
disp(size(var_CG.var_AtA.Dij{1}));
for n_iter = 1:N_iter
    %% Active index update
    [rhs,var_CG.var_AtA] = update_ac(x,var_CG.var_AtA);
    b = AtmX(rhs,var_CG.var_AtA);

    %% Use CG to solve the first subproblem
    [x,cg_err,cg_n] = CG033114(x,b+mup*(zp-up)+mu*Wtx(z-u,wp),var_CG.JtJ,var_CG,var_CG.cg_tol,var_CG.cg_iter);
    disp(norm(x));
    %% Update the second subproblem
    xmin = mu_min;
    xt = up + x;
    id = find(xt >= xmin/2);
    zp = zeros(nX,1,'single');
    zp(id) = max(xt(id),xmin);

    %% Update the first Lagrangian multiplier
    up = up + x - zp;
    
    %% Update the fourth subproblem
    xs = Wx(x,wp);%
    xt = u + xs;%
    [tmp,id] = sort(xt,'descend');
    z = zeros(ip.ns,1,'single'); 
    z(id(1:N)) = xt(id(1:N));

    %% Update the second Lagrangian multiplier
    u = u + xs - z;

    %% Project to MMU constraint
    xp = x;
    xp(xp < eps) = 0;
    xp(intersect(find(xp >= eps),find(xp < mu_min))) = mu_min;
    
    %% Calculate plan parameter
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

    if n_iter == N_iter
        xs = Wx(xp,wp);
        [xp,Id] = shrink_gs(xp,N,wp);
        sg = Id;
    end
end







