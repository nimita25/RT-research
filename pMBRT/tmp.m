% N_dij=1;
% id_obj=cell(N_obj,N_dij);
% for i=1:N_obj
%     id_obj(i,:)=c(c_obj(i));
% end
% 
% Dij = [];
% for i_id=1:nN
%     for i_ctc=1:numel(ctc)
%         Dij = [Dij all_Dij{i_id,i_ctc}];
%     end
% end
% Dij = {Dij};
% [nY,nX]=size(Dij{1});
% 
% AmX=@AmX_v3;
% AtmX=@AtmX_v3;
% var_plot=struct('n_c',n_c,'dmax',px*1.2);
% mu_i=6;mu_x=-6;mu_z=0;wr=0.1; % parameters
% para=struct('Dij',{Dij},'nX',nX,'nY',nY,'N_c',N_c,'n_c',n_c,'c',{c},'isC',uint32(0),...
%     'N_obj',N_obj,'type_obj',type_obj,'w_obj',w_obj,'N_dij',N_dij,...
%     's_obj',s_obj,'n_obj',n_obj,'id_obj',{id_obj},'c_obj',c_obj,'nN',nN,...
%     'mu_i',mu_i,'mu_x',mu_x,'mu_z',mu_z,'wr',wr, ...
%     'all_Dij',{all_Dij},'id',id,'ctc',ctc);
% ip=struct('N_iter',[],'nplot',10000,'var_plot',var_plot,'isC',uint32(0),'mu_min',mu_min);
% 
% x0=ones([nX 1],'single');
% maxAtA=norm(AtmX(AmX(x0,para),para))/norm(x0);
% 
% mup=0.01;
% mu = 1e-5;
% ip.mup=maxAtA*mup;
% ip.mu = mu;
% ip.N_iter_mip=N_iter_mip;
% ip.N_iter=N_iter;
% 
% wp=struct('Dij',{Dij},'intp',intp,'nX',nX,'nY',nY,'tvx',tvx,'tvz',tvz,...
%     'N_dij',N_dij,'isC',uint32(0));
% % wp=struct('Dij',{Dij},'nX',nX,'nY',nY,...
% %     'N_dij',N_dij,'isC',uint32(0));
% wps=struct('isC',uint32(0));
% 
% if mod(method,2)==1
%         var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(0),...
%             'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
%             'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
%             'wps',wps,'mu_xs',[],'Wxr',@BX_v1,'Wtxr',@BtX_v1,'wpr',wp,'mu_xr',1);
%         tic;[x0,y0]=admm_mmu1_mip(ip,var_CG);toc;
%         n = 0;
%         Dij = [];
%         nBeams = zeros(numel(ctc),1);
%         y1 = zeros(numel(id)*numel(ctc),1);
%         for i_id = 1:numel(id)
%             %[~,tmp_id] = max(y0(n+(1:numel(ctc))));
%             if i_id == 1 || i_id == 4
%                 tmp_id = 1;
%             else
%                 tmp_id = 1;
%             end
%             y1(n+tmp_id) = 1;
%             Dij = [Dij var_CG.var_AtA.all_Dij{i_id,tmp_id}];
%             nBeams(i_id) = size(var_CG.var_AtA.all_Dij{i_id,tmp_id},2);
%             n = n+numel(ctc);
%         end
%         disp(y1);
%         Dij = {Dij};
%         [nY,nX]=size(Dij{1});
%         disp([nY,nX]);
%         var_CG.var_AtA.Dij = Dij;
%         var_CG.var_AtA.nY = nY;
%         var_CG.var_AtA.nX = nX;
%         tic;x0=admm_mmu1(ip,var_CG);toc;
% else
%         var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(1),...
%             'JtJ','JtJ_wtw','cg_iter',10,'cg_tol',1e-5,...
%             'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',@wx_id,'Wtxs',@wx_id,...
%             'wps',wps,'mu_xs',[],'Wxr',@BX_v3,'Wtxr',@BtX_v3,'wpr',wp,'mu_xr',1);
%         %tic;[x0,y0]=admm_mmu4_mip(ip,var_CG);toc;
%         n = 0;
%         Dij = [];
%         nBeams = zeros(numel(ctc),1);
%         y1 = zeros(numel(id)*numel(ctc),1);
%         for i_id = 1:numel(id)
%             %[~,tmp_id] = max(y0(n+(1:numel(ctc))));
%             if i_id == 1 || i_id == 4
%                 tmp_id = 3;
%             else
%                 tmp_id = 2;
%             end
%             y1(n+tmp_id) = 1;
%             Dij = [Dij var_CG.var_AtA.all_Dij{i_id,tmp_id}];
%             nBeams(i_id) = size(var_CG.var_AtA.all_Dij{i_id,tmp_id},2);
%             n = n+numel(ctc);
%         end
%         disp(y1);
%         Dij = {Dij};
%         [nY,nX]=size(Dij{1});
%         disp([nY,nX]);
%         var_CG.var_AtA.Dij = Dij;
%         var_CG.var_AtA.nY = nY;
%         var_CG.var_AtA.nX = nX;
%         var_CG.wpr.Dij = Dij;
%         var_CG.wpr.nY = nY;
%         var_CG.wpr.nX = nX;
%         tic;x0=admm_mmu4(ip,var_CG);toc;
% end
% 
% % 4. evaluation
% xp=double(x0);
% xp(xp<eps)=0;
% xp(intersect(find(xp>=eps),find(xp<mu_min)))=mu_min;
% 
% if 1
% n_ctv98=round(n_c(1)*0.95);
% y0=Dij{1}*double(xp);
% y=y0(c{1});
% y2=sort(y,'descend');
% factor=px/y2(n_ctv98);
% x=xp*factor;
% else
% x=xp;
% end
% 
% d=Dij{1}*double(x);
% 
% %% Get output parameters
% % Calculate objective function value
% [obj,D]=calc_obj_dvh(x,var_CG.var_AtA);
% ObjFnVal = mean(sum(obj));
% 
% %Dmax calculation
% y=d(c{1});
% y2=sort(y,'descend');
% Dmax=y2(1)/px; 
% 
% %CI calculation
% V100 = sum(y2>=px);
% V = numel(c{1});
% V100_all = sum(d>=px);
% CI = (V100*V100)/(V*V100_all); 
% 
% %Mean dose in target, body, OAR calculation
% MD = zeros(N_c,1);
% for i_N = 1:N_c
%     y = d(c{i_N});
%     MD(i_N) = mean(y)*100/px;
% end
% 
% % Calculate mean doses in in target, body, OAR per angle
% tmp_index = 0;
% mean_doses = zeros(numel(id),N_c);
% for i = 1:numel(id)
%     dd = Dij{1}(:,tmp_index+1:tmp_index+nBeams(i))*double(x(tmp_index+1:tmp_index+nBeams(i)));
%     tmp_index = tmp_index+nBeams(i);
%     for j = 1:N_c
%         mean_doses(i,j) = mean(dd(c{j}));
%     end
% end
% disp(mean_doses);
% 
% %% Generate 3D image
% % if method>1
% d3d=reshape(d,doseGrid.dimensions);
% % else
% % d=reshape(d,ct.cubeDim);
% % end
% %figure;imshow3D(d3d,[0,px]);
% 
% %% Save output
% outp.x = x;
% outp.d = d;
% outp.d3d = d3d;
% outp.factor = factor;
% outp.Dmax = Dmax;
% outp.CI = CI;
% outp.MD = MD;
% outp.mean_doses = mean_doses;
% outp.ctc = ctc;
% outp.id = id;
% outp.y_bin = reshape(y1,[numel(ctc),numel(id)]);
% outp.ObjFnVal = ObjFnVal;
% 
% for i=1:numel(id)
%     if id(i)==45
%         dk=intp_45*d;
%         D_45 = mean(dk)*100/px;
%         outp.D_45 = D_45;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_45 = PVDR;
%     elseif id(i) == 135
%         dk=intp_135*d;
%         D_135 = mean(dk)*100/px;
%         outp.D_135 = D_135;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_135 = PVDR;
%     elseif id(i) == 240
%         dk=intp_240*d;
%         D_225 = mean(dk)*100/px;
%         outp.D_225 = D_225;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_225 = PVDR;
%     else
%         dk=intp_315*d;
%         D_315 = mean(dk)*100/px;
%         outp.D_315 = D_315;
%         dk1=sort(dk,'descend');
%         dk2=dk1(dk1>0.01);
%         d10 = ceil(0.1*numel(dk2));
%         d80 = ceil(0.8*numel(dk2));
%         PVDR = dk2(d10)/dk2(d80);
%         outp.PVDR_315 = PVDR;
%     end
% end
% 
% if mod(method,2)==1
%     save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(numel(ctc)) '_' num2str(method) '.mat'],"-struct",'outp');
%     %save([ptid '\res_' ptid '_' 'tmp.mat'],"-struct",'outp');
%     %save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '.mat'],'x','d','d3d','factor', 'mean_doses');
% else
%     %save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(method) '_' num2str(mu_i) '_' num2str(-mu_x) '_' num2str(-mu_z) '_' num2str(wr) '.mat'],'x','d','factor', 'mean_doses');
%     outp.mu_i = mu_i;
%     outp.mu_x = mu_x;
%     outp.mu_z = mu_z;
%     outp.wr = wr;
%     save([ptid '\res_' ptid '_' num2str(numel(id)) '_' num2str(numel(ctc)) '_' num2str(method) '.mat'],"-struct",'outp');
% end
% 
% 
px = 2.12;
%load('C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\2286842_4_tvm.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\2286842_4_intp.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_2_ctc3333.mat')
%load('C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_3_ctc4334.mat')

% for i=1:numel(id)
%     if id(i)==45
%         dk=intp_45*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_45 = D_bev;
%         outp.PVDR_45 = [pvdr1,pvdr2];
%     elseif id(i) == 135
%         dk=intp_135*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_135 = D_bev;
%         outp.PVDR_135 = [pvdr1,pvdr2];
%     elseif id(i) == 225
%         dk=intp_225*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_225 = D_bev;
%         outp.PVDR_225 = [pvdr1,pvdr2];
%     else
%         dk=intp_315*d;
%         [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
%         outp.D_315 = D_bev;
%         outp.PVDR_315 = [pvdr1,pvdr2];
%     end
% end


file1 = {'C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_2_ctc3333.mat',
    'C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_2_ctc3443.mat'
    'C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_2_ctc3553.mat'
    'C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_2_ctc4444.mat'
    'C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\res_2286842_4_3_2_ctc5335.mat'};
for j = 1:5
file = file1{j};%'C:\Users\nshinde\Desktop\pMBRT\2286842\1243050\res_1243050_3_3_2_ctc335.mat';
disp(file)
load(file)
load('C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\2286842_4_intp.mat')
load('C:\Users\nshinde\Desktop\pMBRT\2286842\2286842\2286842_4_tvm.mat')
clear outp

outp.x = x;
outp.d = d;
outp.d3d = d3d;
outp.factor = factor;
outp.Dmax = Dmax;
outp.CI = CI;
outp.MD = MD;
outp.mean_doses = mean_doses;
outp.ctc = ctc;
outp.id = id;
outp.y_bin = y_bin;
outp.ObjFnVal = ObjFnVal;
outp.mu_i = mu_i;
outp.mu_x = mu_x;
outp.mu_z = mu_z;
outp.wr = wr;

for i=1:numel(id)
    if id(i)==45
        dk=intp_45*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_45 = D_bev;
        outp.PVDR_45 = [pvdr1,pvdr2];
    elseif id(i) == 135
        dk=intp_135*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_135 = D_bev;
        outp.PVDR_135 = [pvdr1,pvdr2];
    elseif id(i) == 225
        dk=intp_225*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_225 = D_bev;
        outp.PVDR_225 = [pvdr1,pvdr2];
    else
        dk=intp_315*d;
        [D_bev, pvdr1, pvdr2] = calc_bev_para(dk,px);
        outp.D_315 = D_bev;
        outp.PVDR_315 = [pvdr1,pvdr2];
    end
end

save(file,"-struct",'outp');
end