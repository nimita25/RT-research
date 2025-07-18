classdef functionsContainer
   methods

       function [objFn] = define_BED_max_objective(obj,zz,N_obj,tmp_Dij,i,nm,lambd,x0,unique_fields)
        objFn = 0;
        n = 0;
        for nf = 1:N_obj
            tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
            objFn = objFn+(zz(nf)-tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))+lambd{nm}(i,nf))^2;
            n = n+tmp_index;
        end
        end
        
        function [objFn] = define_BED_mean_objective(obj,zz,N_obj,tmp_Dij,nm,lambd,x0,n_c,unique_fields)
        objFn = 0;
        n = 0;
        for nf = 1:N_obj
            tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
            objFn = objFn+sum((zz(:,nf)-tmp_Dij{mod(nf-1,unique_fields)+1,nm}*x0(n+ (1:tmp_index))+lambd{nm}(:,nf)).^2);
            n = n+tmp_index;
        end
        % for i = 1:n_c
        %     n = 0;
        %     for nf = 1:N_obj
        %         tmp_index = size(tmp_Dij{mod(nf-1,unique_fields)+1,nm},2);
        %         objFn = objFn+(zz(i,nf)-tmp_Dij{mod(nf-1,unique_fields)+1,nm}(i,:)*x0(n+ (1:tmp_index))+lambd{nm}(i,nf))^2;
        %         n = n+tmp_index;
        %     end
        % end
        end
        
        function [cineq, ceq] = define_BED_max_constraint(obj,zz,N_obj,nm,rho,BED_max,RBE)
        cineq = 0;
        for nf = 1:N_obj
            cineq = cineq+RBE*(zz(nf))+rho(nm)*(zz(nf))^2;
        end
        cineq = cineq - BED_max;
        %disp(cineq);
        
        ceq=[];
        end
        
        function [cineq, ceq] = define_BED_mean_constraint(obj,zz,N_obj,nm,rho,BED_max,n_c,RBE)
        cineq = 0;
        for nf = 1:N_obj
            cineq = cineq+sum(RBE*(zz(:,nf))+rho(nm)*(zz(:,nf)).^2);
        end
        % for i = 1:n_c
        %     for nf = 1:N_obj
        %         cineq = cineq+RBE*(zz(i,nf))+rho(nm)*(zz(i,nf))^2;
        %     end
        % end
        cineq = cineq - n_c*BED_max(nm);
        
        ceq=[];
        end
        
        function [id_min,Viol_DVHmin] = find_active_indices_DVH_min(obj,x,DVH_min_perc,Cost_matrix,DVH_min)
        Viol_DVHmin = 1;
        [eval2, I] = sort(x,'descend');
        ai = ceil(DVH_min_perc*size(Cost_matrix,1));
        if eval2(ai) < DVH_min
            Viol_DVHmin = min(Viol_DVHmin, eval2(ai)/DVH_min);
            id1 = find(eval2<DVH_min,1,'first');
            if isempty(id1)
                id_min = I(1:ai);
            else
                id_min = I(id1:ai);
            end
        end
        end

        function [id,Viol_BEDDVH] = find_active_indices_DVH_max(obj,DVH_eval,num_DVH,DVH_perc,n_c,BED_DVH_max,id)
        Viol_BEDDVH = 0;
        [eval2, I] = sort(DVH_eval,'descend');
        for n = 1:num_DVH
            ai = ceil(DVH_perc(n)*n_c);
            if eval2(ai) > BED_DVH_max(n)
                Viol_BEDDVH = max(Viol_BEDDVH, eval2(ai)/BED_DVH_max(n));
                id1 = find(eval2<=0.99*BED_DVH_max(n),1,'first');
                if isempty(id1)
                    id{n} = I(ai:end);
                else
                    id{n} = I(ai:id1);
                end
            end
        end
        end

        function [result] = JtJ(obj,J)
            [m,n]=size(J);
            k=n-1;
            kc=k+1;
            %tic;
            B=zeros(n);
            B(:,kc)=sum(J.^2);
            for i=1:k
            tmp=sum(J(:,1:end-i).*J(:,i+1:end));
            B(1:end-i,kc-i)=tmp;
            B(i+1:end,kc+i)=tmp;
            end
            result=spdiags(B,-k:k,n,n);
            %toc;
        end
   end
end

