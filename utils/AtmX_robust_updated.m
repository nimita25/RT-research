function x=AtmX_robust_updated(y,var)
Dij=var.Dij;
id_obj=var.id_obj;
N_obj=var.N_obj;
w_obj=var.w_obj;
n_c=var.n_c;
c_obj=var.c_obj;
nY=var.nY;
nX=var.nX;
N_Dij=var.N_Dij;

%x=zeros([nX,1]);
x = zeros([size(Dij{1},2),1]); %added by Nimita
for j=1:N_Dij
    y0=zeros([nY,1]);
    for i=1:N_obj
        y0(id_obj{i,j})=y0(id_obj{i,j})+w_obj(i,j)/n_c(c_obj(i))*y{i,j};
    end
    %disp(size(y0));
    %new_y0 = []; %added by Nimita
    %for i=1:N_obj  %added by Nimita
    %    new_y0=[new_y0;y0(id_obj{i,j})];  %added by Nimita
    %end  %added by Nimita
    %disp(size(Dij{j}));
    %disp(size(new_y0));
    x=x+Dij{j}'*y0; 
    %x=x+Dij{j}'*new_y0; %added by Nimita
end
x=single(x);
