function [rhs,var]=update_ac_robust(x,var)
Dij=var.Dij;
c=var.c;
N_obj=var.N_obj;
n_obj=var.n_obj;
s_obj=var.s_obj;
c_obj=var.c_obj;
type_obj=var.type_obj;
N_dij=var.N_Dij;

% flag_obj=zeros([N_obj N_dij]);
rhs=cell(N_obj,N_dij);
id_obj=cell(N_obj,N_dij);

for j=1:N_dij
y0=Dij{j}*double(x);

% id_obj=var.id_obj;
% sgn_obj=var.sgn_obj;
for i=1:N_obj
    if type_obj(i)==1 % max constraint
        id=c{c_obj(i)};
        y=y0(id);
        [y2,I]=sort(y,'descend');
        id=id(I);
        if y2(n_obj(i))>s_obj(i)
%             flag_obj(i,j)=1;
            id1=find(y2<=1.01*y2(n_obj(i)),1,'first');
            id2=find(y2<=0.99*s_obj(i),1,'first');
            id_obj{i,j}=id(id1:id2);                
%             id2=find(y2<=s_obj(i),1,'first');
%             id_obj{i}=id(n_obj(i):id2);
            rhs{i,j}=ones([numel(id_obj{i,j}) 1])*s_obj(i);
%             disp(['C' num2str(j) '-' num2str(i) 'dvhmax ' num2str([s_obj(i) y2(n_obj(i)) numel(id_obj{i,j})/numel(id)])])
        end
    elseif type_obj(i)==2 % min constraint
        id=c{c_obj(i)};
        y=y0(id);
        [y2,I]=sort(y,'descend');
        id=id(I);
%         disp(num2str([y2(n_obj(i)) s_obj(i)]))
        if y2(n_obj(i))<s_obj(i)
%             flag_obj(i,j)=1;
%             id1=find(y2<s_obj(i),1,'first');
%             id2=find(y2<=0.99*y2(n_obj(i)),1,'first');
%             id_obj{i,j}=id(id1:id2);
            id1=find(y2<s_obj(i),1,'first');
            id_obj{i,j}=id(id1:n_obj(i));
            rhs{i,j}=ones([numel(id_obj{i,j}) 1])*s_obj(i);
%             disp(['C' num2str(j) '-' num2str(i) 'dvhmin ' num2str([s_obj(i) y2(n_obj(i)) numel(id_obj{i,j})/numel(id)])])
        end
    elseif type_obj(i)==3 % max ptv/body
        id=c{c_obj(i)};
        y=y0(id);
        id2=find(y>s_obj(i));
        if ~isempty(id2)
%             flag_obj(i,j)=1;
            id_obj{i,j}=id(id2);
            rhs{i,j}=ones([numel(id_obj{i,j}) 1])*s_obj(i);
%             disp(['C' num2str(j) '-' num2str(i) 'max ' num2str([s_obj(i) max(y(id2)) mean(y(id2)) numel(id_obj{i,j})/numel(id)])])
        end
    elseif type_obj(i)==0
%         flag_obj(i,j)=1;
        id_obj{i,j}=c{c_obj(i)};
        rhs{i,j}=ones([numel(id_obj{i,j}) 1])*s_obj(i);
%         disp(['C' num2str(i) 'mean ' num2str([s_obj(i)])])
    end
end
end


% %added by Nimita
% Dij_new = [];
% for j=1:N_dij
%     for i=1:N_obj
%         Dij_new = [Dij_new;Dij{1}(id_obj{i,j},:)];
%     end
% end
% var.Dij{1} = Dij_new;


var.id_obj=id_obj;
% var.flag_obj=flag_obj;

