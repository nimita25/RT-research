function x = admm_mmu_2601(ip, var_CG)
% AtAmX Calculates the result of A' * A * X.
%
% INPUT:
%   ip      - Structure containing relevant information for plot.
%   var_CG  - Structure containing relevant information for model and CG.
%
% OUTPUT:
%   x - optimized beam intensity

% Extract parameters from input
AtmX = var_CG.AtmX; % A^T
AmX = var_CG.AmX; % A
N_dij = numel(var_CG.var_AtA.Dij); % number of cases, set 1 as default
n_ctv95 = round(var_CG.var_AtA.n_c(1) * 0.95);
nX = var_CG.var_AtA.nX;
px = var_CG.var_AtA.s_obj(1); % prescription dose for the target

%========================================================
mu_min = ip.mu_min; % MMU threshold
N_iter = ip.N_iter; % Number of iterations
mup = ip.mup; % Augmented Lagrangian Parameter
var_CG.mu_xs = ip.mup; % parameter for the first subproblem
nplot = ip.nplot;
var_plot = ip.var_plot;

% Initialize variables
x = zeros([nX 1], 'single'); % primal variable
zp = x; % primal variable
up = x; % dual variable

% ICR and ADMM iterations
for n_iter = 1:N_iter

    % Update DVH active index
    [rhs, var_CG.var_AtA] = update_ac(x, var_CG.var_AtA); %%
    b = AtmX(rhs, var_CG.var_AtA);

    % 1. Update x via L2 problem using CG method
    [x, cg_err, cg_n] = CG033114(x, b + mup * (zp - up), var_CG.JtJ, var_CG, var_CG.cg_tol, var_CG.cg_iter);

    % 2. Update second primal variable zp
    xt = up + x;
    id = find(xt >= mu_min / 2);
    zp = zeros(nX, 1, 'single');
    zp(id) = max(xt(id), mu_min);

    % 3. Update multipliers
    up = up + x - zp;

    % Project x to MMU constraint
    xp = x;
    xp(xp < mu_min / 2) = 0;
    xp(intersect(find(xp >= mu_min / 2), find(xp < mu_min))) = mu_min;

    % Calculate objective and DVH metrics
    [obj, D] = calc_obj_dvh(xp, var_CG.var_AtA);
    obj_total = sum(obj);

    D95 = zeros(N_dij, 1);
    Dmax = zeros(N_dij, 1);
    for i = 1:N_dij
        y2 = sort(D{1, i}, 'descend');
        D95(i) = y2(n_ctv95) / px;
        Dmax(i) = y2(1) / px;
    end

    % Plot DVH if needed
    if mod(n_iter, nplot) == 0
        plotdvh(D, var_plot);
    end

    % Display iteration results
    disp(num2str([n_iter mean(obj_total) mean(D95) mean(Dmax)]))
end
end