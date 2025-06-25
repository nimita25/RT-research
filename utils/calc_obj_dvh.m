function [obj, d] = calc_obj_dvh(x, var)
% calc_obj_dvh calculates the objective values.
%
% INPUT:
%   x   - Beam intensity vector.
%   var - Structure containing relevant information.
%
% OUTPUT:
%   obj - Objective values for each object.
%   d   - Dose distributions for different structure

% Extract necessary information from the input structure
Dij = var.Dij{1};
N_obj = var.N_obj;
id_obj = var.id_obj;
w_obj = var.w_obj;
c_obj = var.c_obj;
s_obj = var.s_obj;
N_c = var.N_c;
n_c = var.n_c;
c = var.c;

% calculate dose vector
y0 = Dij * double(x);

% Calculate objective values for each object
obj = zeros(N_obj, 1);
for i = 1:N_obj
    obj(i) = w_obj(i) / n_c(c_obj(i)) * sum((y0(id_obj{i}) - s_obj(i)).^2);
end

% Extract dose distributions for different structure
d = cell(N_c, 1);
for i = 1:N_c
    d{i} = y0(c{i});
end
end
