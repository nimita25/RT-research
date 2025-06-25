function [x,y]=admm_mmu1_mip(ip,var_CG)

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

x=zeros([nX 1],'single');up=x;zp=up;

y=ones(numel(id)*numel(ctc),1)/numel(ctc);
%y=zeros(numel(id)*numel(ctc),1);
%y(:,2) =1;
%disp(y);
mu = ip.mu;
u = ones(numel(id)*numel(ctc),1);
N_obj = var_CG.var_AtA.N_obj;

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

for n_iter=1:N_iter_mip


    

    [rhs,var_CG.var_AtA]=update_ac(x,var_CG.var_AtA); %%
    b=AtmX(rhs,var_CG.var_AtA);


    % 1. update x via L2 problem %
    [x,cg_err,cg_n]=CG033114(x,b+mup*(zp-up),var_CG.JtJ,var_CG,var_CG.cg_tol,var_CG.cg_iter);

    % 2. mmu constraint %
    xmin=mu_min;
    xt=up+x;
    i_d=find(xt>=xmin/2);
    zp=zeros(nX,1,'single');
    zp(i_d)=max(xt(i_d),xmin);
    up=up+x-zp;

    % 3. update binary variable %
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

    [rhs,var_CG.var_AtA]=update_ac(y,var_CG.var_AtA); %%
    b=AtmX(rhs,var_CG.var_AtA);
    ATA = [];
    for i_id=1:numel(id)
        for i_ctc=1:numel(ctc)
            tmp_id = (i_id-1)*numel(ctc)+i_ctc;
            tmp_rhs = cell(N_obj,N_dij);
            for j = 1:N_dij
                for i = 1:N_obj
                    tmp_rhs{i,j} = Ax{1}(var_CG.var_AtA.id_obj{i,j},tmp_id);
                end
            end
            %tmp_rhs{1,1} = Ax{1}(:,tmp_id);
            ATA=[ATA AtmX(tmp_rhs,var_CG.var_AtA)];
        end
    end
    BLKD = ones(numel(ctc),numel(ctc));
    for i_id = 1:numel(id)-1
        BLKD = blkdiag(BLKD,ones(numel(ctc),numel(ctc)));
    end
    A_lhs = ATA+mu*BLKD;
    b_rhs = b+mu*(ones(numel(ctc)*numel(id),1)-u);
    y = mldivide(A_lhs,b_rhs);
    %disp(y);

    % 4. Update dual variable corresponding to y %
    n = 0;
    for i_id=1:numel(id)
        u(n+(1:numel(ctc))) = sum(y(n+(1:numel(ctc))))-1;
        n = n+numel(ctc);
    end

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