function [xp,sg,s_large] = admm_ip(ip,var_CG,NNZ,x_init)
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
Wtx = var_CG.Wtxr;
Wx = var_CG.Wxr;
wp = var_CG.wpr;
x = x_init;
up = x; zp = up;
%Initialize s
s = Wx(x,wp);
[tmp,id] = sort(s,'descend');
s = zeros(ip.ns,1,'single'); 
s(id(1:N)) = 1;
u = s; z = u;
sg = cell(N,1);
n_iter = 0;
N_iter = ip.N_iter;

%Save original Dij and use Dij.xj and Dij.sj, where . is Hadamard product
Dij_orig = var_CG.var_AtA.Dij;
size_Dij = size(Dij_orig{1});

%% Optimization
for n_iter = 1:N_iter
    %% Active index update
    % Update Dij.sj
    var_CG.var_AtA.Dij = Dij_orig;
    [rhs,var_CG.var_AtA] = update_ac(x,var_CG.var_AtA); %get active columns of Dij
    s_large = Wtx(z,wp);%Wtx(s,wp); %This is the information about active energy layers
    %with next 2 lines, only the columns of Dij corresponding to active energy layers are used
    %This is defined specifically for D (times) s
    %tmp_AtA = var_CG.var_AtA; %create a temp copy
    var_CG.var_AtA.Dij{1} = double(s_large)'.*var_CG.var_AtA.Dij{1}; %temp info about D: D (times) s
    b = AtmX(rhs,var_CG.var_AtA); %Get ((D times s)^T.W)b
    %disp(size(var_CG.var_AtA.Dij{1}))
    
    %% Update x - Solve linear system of equations
    Dij_new = [];
    for j=1:N_dij
        for i=1:N_obj
            Dij_new = [Dij_new;var_CG.var_AtA.Dij{1}(var_CG.var_AtA.id_obj{i,j},:)];
        end
    end
    %disp(size(Dij_new))
    var_CG.var_AtA.Dij{1} = Dij_new;
    A = var_CG.var_AtA.Dij{1}'*var_CG.var_AtA.Dij{1}+double(mup)*speye(size_Dij(2));%probably need to update this var_CG.Dij
    b_lin = b+mup*(zp-up);
    b_lin = sparse(double(b_lin));
    x = mldivide(A,b_lin);
    x = full(x);
    disp(norm(x));

    %% Update z1 - the second subproblem
    xmin = mu_min;
    xt = up + x;
    id = find(xt >= xmin/2);
    zp = zeros(nX,1,'single');
    zp(id) = max(xt(id),xmin);

    %% Update u1 - the first Lagrangian multiplier
    up = up + x - zp;

    %% Update s - the binary variable
    var_CG.var_AtA.Dij = Dij_orig;
    [rhs,var_CG.var_AtA] = update_ac(x,var_CG.var_AtA); %get active columns of Dij
    %with next 2 lines, we define D (times) x
    var_CG.var_AtA.Dij{1} = double(x)'.*var_CG.var_AtA.Dij{1}; %temp info about D: D (times) x
    b = AtmX(rhs,var_CG.var_AtA); %Get ((D times x)^T.W)b

    %% Update s - Solve linear system of equations
    Dij_new = [];
    for j=1:N_dij
        for i=1:N_obj
            Dij_new = [Dij_new;var_CG.var_AtA.Dij{1}(var_CG.var_AtA.id_obj{i,j},:)];
        end
    end
    var_CG.var_AtA.Dij{1} = Dij_new;
    A = var_CG.var_AtA.Dij{1}'*var_CG.var_AtA.Dij{1}+double(mu)*speye(size_Dij(2));%probably need to update this var_CG.Dij
    b_lin = b+mu*(Wtx(z,wp)-Wtx(u,wp));
    b_lin = sparse(double(b_lin));
    A = full(A);
    A_reduced = zeros(length(b_lin),ip.ns);
    for i = 1:length(b_lin)
        A_reduced(i,:) = Wx(A(i,:)',wp);
    end
    s = mldivide(A_reduced,b_lin);

    %% Update fourth variable
    st = u + s;%
    [tmp,id] = sort(st,'descend');
    z = zeros(ip.ns,1,'single'); 
    z(id(1:N)) = st(id(1:N));

    %% Update the second Lagrangian multiplier
    u = u + s - z;

    %% Update z2 - the fourth subproblem
    % st = u + s;%
    % [tmp,id] = sort(st,'descend');
    % z = zeros(ip.ns,1,'single'); 
    % for idx = id(1:N)
    %     if st(idx) > 0.5
    %         z(idx) = 1;
    %     end
    % end
    % %z(id(1:N)) = st(id(1:N));%change this so that z has integral values
    
    %% Project to MMU constraint
    xp = x;
    xp(xp < eps) = 0;
    xp(intersect(find(xp >= eps),find(xp < mu_min))) = mu_min;
    
    %% Calculate plan parameter
    var_CG.var_AtA.Dij = Dij_orig;
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







