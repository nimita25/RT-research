function x=admm_mmu1(ip,var_CG)

mu_min=ip.mu_min;
eps=mu_min/2;
N_dij=numel(var_CG.var_AtA.Dij);
n_ctv98=round(var_CG.var_AtA.n_c(1)*0.95);
nX=var_CG.var_AtA.nX;
px=var_CG.var_AtA.s_obj(1);
mup=ip.mup;
var_CG.mu_xs=ip.mup;
N_iter=ip.N_iter;AtmX=var_CG.AtmX;AmX=var_CG.AmX;
nplot=ip.nplot;var_plot=ip.var_plot;

x=zeros([nX 1],'single');up=x;zp=up;

for n_iter=1:N_iter

    [rhs,var_CG.var_AtA]=update_ac(x,var_CG.var_AtA); %%
    b=AtmX(rhs,var_CG.var_AtA);
    
    
    % 1. update x via L2 problem %
    [x,cg_err,cg_n]=CG033114(x,b+mup*(zp-up),var_CG.JtJ,var_CG,var_CG.cg_tol,var_CG.cg_iter);
    
    % 2. mmu constraint %
    xmin=mu_min;
    xt=up+x;
    id=find(xt>=xmin/2);
    zp=zeros(nX,1,'single');
    zp(id)=max(xt(id),xmin);
    up=up+x-zp;

    % calc plan parameters %
    xp=x;
    xp(xp<eps)=0;        
    xp(intersect(find(xp>=eps),find(xp<mu_min)))=mu_min;        

    [obj,D]=calc_obj_dvh(xp,var_CG.var_AtA);
    obj_total=sum(obj);
    
    D98=zeros(N_dij,1);
    Dmax=zeros(N_dij,1);
    for i=1:1
    y2=sort(D{1,i},'descend');
    D98(i)=y2(n_ctv98)/px;
    Dmax(i)=y2(1)/px;
    end

    if mod(n_iter,nplot)==0
    plotdvh(D,var_plot);
    end
    disp(num2str([n_iter mean(obj_total) mean(D98) mean(Dmax)]))
end