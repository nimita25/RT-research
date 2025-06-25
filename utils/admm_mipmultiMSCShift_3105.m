function [xp,y] = admm_mipmultiMSCShift_3105(ip,var_CG,x_init,N,qc)
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
var_CG.mu_xr = ip.mu;
wp = var_CG.wpr;
n_gs = wp.nx;
n_ba = single(wp.na);
x = x_init;
up = x; zp = up;
%Initialize y
%y = ones(n_ba,1,'single')/n_ba;
y = ones(n_ba,1,'single')*N/n_ba;%ones(n_ba,1,'single')*N/n_ba;
%y = zeros(n_ba,1);y(3)=1;y(7)=1;y(11)=1;y(15)=1;
u = sum(y) - N;
%u = y;
% z = zeros(n_ba,1,'single');
% if nargin == 5
%     z(isg) = 1;
% else
%     z(1:N) = 1;
% end
% u = z; 
% sg = cell(N,1);
N_iter = ip.N_iter;

%Save original Dij and use Dij.xj and Dij.yj, where . is Hadamard product
Dij_orig = var_CG.var_AtA.Dij;
size_Dij = size(Dij_orig{1});

%% Optimization 
for n_iter = 1:N_iter
    %% Update Dij according to active shifts
    nn = 0;
    for ii = 1:n_ba
        var_CG.var_AtA.Dij{1}(:,nn+(1:n_gs(ii))) = Dij_orig{1}(:,nn+(1:n_gs(ii)))*double(y(ii));
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
    i_d = find(xt >= xmin/2);
    zp = zeros(nX,1,'single');
    zp(i_d) = max(xt(i_d),xmin);

    %% Calculate D(times)x
    Dx = [];
    nn = 0;
    for ii = 1:n_ba
        Dx = [Dx Dij_orig{1}(:,nn+(1:n_gs(ii)))*double(x(nn+(1:n_gs(ii))))];
        nn = nn+n_gs(ii);
    end
    Dx = {Dx};
    var_CG.var_AtA.Dij = Dx;
    [nY, nX] = size(Dx{1});
    var_CG.var_AtA.nX = nX;
    var_CG.var_AtA.nY = nY;

    %% Calculate quantities needed to update s
    %[rhs,var_CG.var_AtA]=update_ac(y,var_CG.var_AtA); %%
    b=AtmX(rhs,var_CG.var_AtA);
    DTD = [];
    for ii = 1:n_ba
        tmp_rhs = cell(N_obj,N_dij);
        for j = 1:N_dij
            for i = 1:N_obj
                tmp_rhs{i,j} = Dx{1}(var_CG.var_AtA.id_obj{i,j},ii);
            end
        end
        DTD=[DTD AtmX(tmp_rhs,var_CG.var_AtA)];
    end
       
    
    
    %% Update y -> the variable that chooses active shifts
    A_lhs = DTD+mu*ones(n_ba,n_ba);
    b_rhs = b+mu*(N-u);
    disp(b');
    disp(b_rhs');
    % disp(DTD);
    % disp(mu);
    % disp(A_lhs);
    % disp(b_rhs);
    if qc == 0
        y = mldivide(A_lhs,b_rhs);
    elseif qc == 1 || qc == 3
        b_rhs = -2*b_rhs;
        qprob = qubo(double(A_lhs),double(b_rhs));

        tic;
        % opts = optimset(MaxIter=500);
        % qa = qaoa(NumLayers=1,NumShots=100,OptimizationSolverOptions=opts);
        % result = solve(qprob,Algorithm=qa);
        ts=tabuSearch(MaxStallTime=10);
        result = solve(qprob,Algorithm=ts);
        % result = solve(qprob);
        toc

        y = result.BestX;
        %y(7)=1;y(8)=0;

    end
    disp(y');



    %% Update the Lagrangian multipliers
    up = up + x - zp;
    u = u + sum(y) - N;
    disp(u');

    %% Project to MMU constraint
    xp = x;
    xp(xp < eps) = 0;
    xp(intersect(find(xp >= eps),find(xp < mu_min))) = mu_min;

    %% Update Dij according to active shifts
    nn = 0;
    Dx = [];
    for ii = 1:n_ba
        Dx = [Dx Dij_orig{1}(:,nn+(1:n_gs(ii)))*double(y(ii))];
        nn = nn+n_gs(ii);
    end
    Dx = {Dx};
    var_CG.var_AtA.Dij = Dx;
    [nY, nX] = size(Dx{1});
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


end





