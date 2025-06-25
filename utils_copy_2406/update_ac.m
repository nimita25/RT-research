function [rhs, var] = update_ac(x, var)
% update_ac Updates the active index.
%
% INPUT:
%   x   - Beam intensity vector.
%   var - Structure containing relevant information.
%
% OUTPUT:
%   rhs - Cell array of prescription dose for each DVH objective.
%   var - Updated structure.


%% 
Dij = var.Dij{1};
c = var.c;
N_obj = var.N_obj;
n_obj = var.n_obj;
s_obj = var.s_obj;
c_obj = var.c_obj;
type_obj = var.type_obj;
y0 = Dij * double(x);
id_obj = cell(N_obj, 1);
rhs = cell(N_obj, 1);

%% Iterate through each constraint
for i = 1:N_obj
    if type_obj(i) == 1 % max constraint
        id = c{c_obj(i)};
        y = y0(id);
        [y2, I] = sort(y, 'descend');
        id = id(I);
        
        if y2(n_obj(i)) > s_obj(i)
            id1 = find(y2 <= 1.01 * y2(n_obj(i)), 1, 'first');
            id2 = find(y2 <= 0.99 * s_obj(i), 1, 'first');
            id_obj{i} = id(id1:id2);
            rhs{i} = ones([numel(id_obj{i}) 1]) * s_obj(i);
        end
    elseif type_obj(i) == 2 % min constraint
        id = c{c_obj(i)};
        y = y0(id);
        [y2, I] = sort(y, 'descend');
        id = id(I);
        
        if y2(n_obj(i)) < s_obj(i)
            id1 = find(y2 < s_obj(i), 1, 'first');
            id_obj{i} = id(id1:n_obj(i));
            rhs{i} = ones([numel(id_obj{i}) 1]) * s_obj(i);
        end
    elseif type_obj(i) == 3 % max ptv/body
        id = c{c_obj(i)};
        y = y0(id);
        id2 = find(y > s_obj(i));
        
        if ~isempty(id2)
            id_obj{i} = id(id2);
            rhs{i} = ones([numel(id_obj{i}) 1]) * s_obj(i);
        end
    elseif type_obj(i) == 0
        id_obj{i} = c{c_obj(i)};
        rhs{i} = ones([numel(id_obj{i}) 1]) * s_obj(i);
    end
end
var.id_obj = id_obj;
end
