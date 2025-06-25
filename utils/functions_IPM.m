function [obj] = define_obj(u,N_obj,Cost_matrix,Cost_b)
% Define the objective function of BED error minimization problem
% This is the function used by fmincon - the interior point method
%disp('Evaluating objective!')

obj = 0;
n = 0;
for i = 1:N_obj
    tmp_index = size(Cost_matrix{i},2);
    %disp(max(Cost_matrix{i}*u(n+ (1:tmp_index))));
    %disp(min(Cost_matrix{i}*u(n+ (1:tmp_index))));
    obj = obj+norm(Cost_matrix{i}*u(n+ (1:tmp_index))-Cost_b)^2;
    %disp(obj);
    n = n+tmp_index;
end

%Add log penality for dose exceeding 1.1*px in tumor
lambd = 10;
n = 0;
for i = 1:N_obj
    tmp_index = size(Cost_matrix{i},2);
    obj = obj-lambd*sum(log(1.1*Cost_b-Cost_matrix{i}*u(n+ (1:tmp_index))));
    %disp(obj);
    n = n+tmp_index;
end
%obj = obj-log(sum((u).^2));
end

function [cineq, ceq] = define_BED_max(u,rho,c,tmp_Dij,BED_max,n_c,N_obj,type_constr)
%This function defines the BED max constraint
% This is the function used by fmincon - the interior point method
%disp('Evaluating nonlinear constraints!')
n = 0;
len_cineq = 0;
for i = 1:length(n_c)
    if type_constr(i) == 1
        len_cineq = len_cineq+n_c(i);
    elseif type_constr(i) == 2
        len_cineq = len_cineq+1;
    end
end
cineq = zeros([len_cineq,1]);
for m_c = 3:numel(c)
    if type_constr(m_c-2) == 1
        size_c = size(c{m_c},1);
        m = 0;
        for i = 1:N_obj
            %size_c = size(tmp_Dij{i,m_c-2},1);
            %disp(size(u))
            %disp(m)
            mm = size(tmp_Dij{i,m_c-2},2);
            %disp(mm)
            %disp(size_c)
            %disp(n+size_c)
            cineq(n+ (1:size_c)) = cineq(n+ (1:size_c))+tmp_Dij{i,m_c-2}*u(m+(1:mm)) + rho(m_c-2)*(tmp_Dij{i,m_c-2}*u(m+(1:mm))).^2 - BED_max(m_c-2)*ones(size_c,1)/N_obj;
            %disp(nnz(cineq));
            m = m+mm;
        %cineq = u.^2-2;
        end
        n=n+size_c;
    elseif type_constr(m_c-2) == 2
        size_c = 1;
        m = 0;
        for i = 1:N_obj
            %size_c = size(tmp_Dij{i,m_c-2},1);
            %disp(size(u))
            %disp(m)
            mm = size(tmp_Dij{i,m_c-2},2);
            %disp(mm)
            %disp(size_c)
            %disp(n+size_c)
            cineq(n+ (1:size_c)) = cineq(n+ (1:size_c))+sum(tmp_Dij{i,m_c-2}*u(m+(1:mm))) + rho(m_c-2)*sum((tmp_Dij{i,m_c-2}*u(m+(1:mm))).^2) - size(c{m_c},1)*BED_max(m_c-2)*ones(size_c,1)/N_obj;
            %disp(nnz(cineq));
            m = m+mm;
        %cineq = u.^2-2;
        end
        n=n+size_c;
    end
end
ceq = [];
end



%obj_handle = @(x0)define_obj(x0,N_obj,Cost_matrix,Cost_b);
%disp(obj_handle(x0))

%cineq_handle = (x0)@define_BED_max(x0,rho,c,tmp_Dij,BED_max,n_c,N_obj,type_constr);
%disp(size(cineq_handle(x0)))


% % Define empty arrays for linear constraints, since our problem does not have any
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% ceq = [];
% ub = [];
% lb = zeros([total_beamlets,1]);
% 
% options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',600000);
% [x,fval,exitflag,output] = fmincon(obj_handle,x0,A,b,Aeq,beq,lb,ub,cineq_handle,options);
% disp(fval);
% disp(exitflag);
% disp(output);
