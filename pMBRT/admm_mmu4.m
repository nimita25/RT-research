function x=admm_mmu4(ip,var_CG)

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

mu_i=var_CG.var_AtA.mu_i;
mu_x=var_CG.var_AtA.mu_x;
mu_z=var_CG.var_AtA.mu_z;
wr=var_CG.var_AtA.wr;
wr1=[wr,-wr,-wr];
r_i=mu_i/wr1(1);
r_x=mu_x/wr1(2);
r_z=mu_z/wr1(3);
var_CG.wpr.r_i=r_i;
var_CG.wpr.r_x=r_x;
var_CG.wpr.r_z=r_z;
intp=var_CG.wpr.intp;
tvx=var_CG.wpr.tvx;
tvz=var_CG.wpr.tvz;
BtX=var_CG.Wtxr;
BX=var_CG.Wxr;

ur1=cell(3,N_dij);
for j=1:N_dij
    ur1{1,j}=zeros(size(intp,1),1);
    ur1{2,j}=zeros(size(tvx,1),1);
    ur1{3,j}=zeros(size(tvz,1),1);
end
zr1=ur1;
b_tv=ur1;

x=zeros([nX 1],'single');up=x;zp=up;

for n_iter=1:N_iter

    [rhs,var_CG.var_AtA]=update_ac(x,var_CG.var_AtA);
    b=AtmX(rhs,var_CG.var_AtA);

    for i=1:3
        for j=1:N_dij
            b_tv{i,j}=zr1{i,j}-ur1{i,j};
        end
    end
    btv=BtX(b_tv,var_CG.wpr);
    
    % 1. update x via L2 problem %
    [x,cg_err,cg_n]=CG033114(x,b+mup*(zp-up)+btv,var_CG.JtJ,var_CG,var_CG.cg_tol,var_CG.cg_iter);
    
    % 2. mmu constraint %
    xmin=mu_min;
    xt=up+x;
    id=find(xt>=xmin/2);
    zp=zeros(nX,1,'single');
    zp(id)=max(xt(id),xmin);
    up=up+x-zp;
    
    z_tem=BX(x,var_CG.wpr);
    for i=1:3
        for j=1:N_dij
            b_tv{i,j}=zr1{i,j}-ur1{i,j};
            zr1{i,j}=sign(z_tem{i,j}+ur1{i,j}).*max(abs(z_tem{i,j}+ur1{i,j})-0.5*wr1(i),0);
            ur1{i,j}=ur1{i,j}+z_tem{i,j}-zr1{i,j};
        end
    end

    % calc plan parameters %
    xp=x;
    xp(xp<eps)=0;        
    xp(intersect(find(xp>=eps),find(xp<mu_min)))=mu_min;        

    [obj,D]=calc_obj_dvh(xp,var_CG.var_AtA);
    obj_total=sum(obj);

    D98=zeros(1,1);
    Dmax=zeros(1,1);
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