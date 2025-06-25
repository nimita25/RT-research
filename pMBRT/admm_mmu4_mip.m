function [x,y]=admm_mmu4_mip(ip,var_CG)

mu_min=ip.mu_min;
eps=mu_min/2;
N_dij=numel(var_CG.var_AtA.Dij);
n_ctv98=round(var_CG.var_AtA.n_c(1)*0.98);
nX=var_CG.var_AtA.nX;
px=var_CG.var_AtA.s_obj(1);
mup=ip.mup;
var_CG.mu_xs=ip.mup;
N_iter_mip=ip.N_iter_mip;AtmX=var_CG.AtmX;AmX=var_CG.AmX;
nplot=ip.nplot;var_plot=ip.var_plot;
id = var_CG.var_AtA.id;
ctc = var_CG.var_AtA.ctc;


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

y=ones(numel(id)*numel(ctc),1)/numel(ctc);
%y=zeros(numel(id)*numel(ctc),1);
%y(:,2) =1;
%disp(y);
mu = ip.mu;
u = ones(numel(id)*numel(ctc),1);
N_obj = var_CG.var_AtA.N_obj;

% Update Dij in var_AtA/para -> used to update the variable in fn f(d)
Dij = [];
for i_id=1:numel(id)
    for i_ctc=1:numel(ctc)
        tmp_id = (i_id-1)*numel(ctc)+i_ctc;
        Dij = [Dij var_CG.var_AtA.all_Dij{i_id,i_ctc}*y(tmp_id)];
    end
end
Dij = {Dij};
[nY,nX]=size(Dij{1});
var_CG.var_AtA.Dij = Dij;
var_CG.var_AtA.nY = nY;
var_CG.var_AtA.nX = nX;

% Update Dij in wpr -> used to update the variable in fn g(dk)
% Dij = [];
% for i_id=1:numel(id)
%     for i_ctc=1:numel(ctc)
%         tmp_id = (i_id-1)*numel(ctc)+i_ctc;
%         Dij = [Dij var_CG.wpr.all_Dij{i_id,i_ctc}*y(tmp_id)];
%     end
% end
% Dij = {Dij};
% [nY,nX]=size(Dij{1});
var_CG.wpr.Dij = Dij;
var_CG.wpr.nY = nY;
var_CG.wpr.nX = nX;

% Start the iterations
for n_iter=1:N_iter_mip

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
    i_d=find(xt>=xmin/2);
    zp=zeros(nX,1,'single');
    zp(i_d)=max(xt(i_d),xmin);
    up=up+x-zp;
    
    z_tem=BX(x,var_CG.wpr);
    for i=1:3
        for j=1:N_dij
            b_tv{i,j}=zr1{i,j}-ur1{i,j};
            zr1{i,j}=sign(z_tem{i,j}+ur1{i,j}).*max(abs(z_tem{i,j}+ur1{i,j})-0.5*wr1(i),0);
            ur1{i,j}=ur1{i,j}+z_tem{i,j}-zr1{i,j};
        end
    end

    % update binary variable %
    %Ax = cell(numel(id),numel(ctc));
    Ax = [];
    n = 0;
    for i_id=1:numel(id)
        for i_ctc=1:numel(ctc)
            %disp(i_id)
            tmp_index = size(var_CG.var_AtA.all_Dij{i_id,i_ctc},2);
            %Ax(i_id,i_ctc) = var_CG.var_AtA.all_Dij{i_id,i_ctc}*x(n+(1:tmp_index));
            Ax = [Ax var_CG.var_AtA.all_Dij{i_id,i_ctc}*double(x(n+(1:tmp_index)))];
            n=n+tmp_index;
        end
    end
    Ax = {Ax};
    [nY,nX]=size(Ax{1});
    var_CG.var_AtA.Dij = Ax;
    var_CG.var_AtA.nY = nY;
    var_CG.var_AtA.nX = nX;
    var_CG.wpr.Dij = Ax;
    var_CG.wpr.nY = nY;
    var_CG.wpr.nX = nX;

    [rhs,var_CG.var_AtA]=update_ac(y,var_CG.var_AtA); %%
    b=AtmX(rhs,var_CG.var_AtA);
    for i=1:3
        for j=1:N_dij
            b_tv{i,j}=zr1{i,j}-ur1{i,j};
        end
    end
    btv=BtX(b_tv,var_CG.wpr);
    ATA = [];
    BTB = [];
    for i_id=1:numel(id)
        for i_ctc=1:numel(ctc)
            tmp_id = (i_id-1)*numel(ctc)+i_ctc;
            tmp_rhs = cell(N_obj,N_dij);
            tmp_rhs_B = cell(3,N_dij);
            for j = 1:N_dij
                for i = 1:N_obj
                    tmp_rhs{i,j} = Ax{1}(var_CG.var_AtA.id_obj{i,j},tmp_id);
                end
            end
            for i = 1:N_dij
                tmp_rhs_B{1,i} = intp*Ax{1}(:,tmp_id);
                tmp_rhs_B{2,i} = tvx*tmp_rhs_B{1,i};
                tmp_rhs_B{3,i} = tvz*tmp_rhs_B{1,i};
            end
            BTB = [BTB BtX(tmp_rhs_B,var_CG.wpr)];
            %tmp_rhs{1,1} = Ax{1}(:,tmp_id);
            ATA=[ATA AtmX(tmp_rhs,var_CG.var_AtA)];
        end
    end
    BLKD = ones(numel(ctc),numel(ctc));
    for i_id = 1:numel(id)-1
        BLKD = blkdiag(BLKD,ones(numel(ctc),numel(ctc)));
    end
    A_lhs = ATA+BTB+mu*BLKD;
    b_rhs = b+btv+mu*(ones(numel(ctc)*numel(id),1)-u);
    y = mldivide(A_lhs,b_rhs);
    %disp(y);

    % Update dual variable corresponding to binary variable %
    n = 0;
    for i_id=1:numel(id)
        u(n+(1:numel(ctc))) = sum(y(n+(1:numel(ctc))))-1;
        n = n+numel(ctc);
    end

    % Update Dij in var_AtA/para
    Dij = [];
    for i_id=1:numel(id)
        for i_ctc=1:numel(ctc)
            tmp_id = (i_id-1)*numel(ctc)+i_ctc;
            Dij = [Dij var_CG.var_AtA.all_Dij{i_id,i_ctc}*double(y(tmp_id))];
        end
    end
    Dij = {Dij};
    [nY,nX]=size(Dij{1});
    var_CG.var_AtA.Dij = Dij;
    var_CG.var_AtA.nY = nY;
    var_CG.var_AtA.nX = nX;

    % Update Dij in wpr
    var_CG.wpr.Dij = Dij;
    var_CG.wpr.nY = nY;
    var_CG.wpr.nX = nX;

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