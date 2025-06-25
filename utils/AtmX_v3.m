function x = AtmX_v3(y, var)
% AtmX_v3 Calculates the result of A' * X.
%
% INPUT:
%   y   - Result vector.
%   var - Structure containing information related to A'.
%
% OUTPUT:
%   x - Result of the product of A' and y.

% Extract necessary information from the input structure
Dij = var.Dij{1};
id_obj = var.id_obj;
N_obj = var.N_obj;
w_obj = var.w_obj;
n_c = var.n_c;
c_obj = var.c_obj;
nY = var.nY;

y0 = zeros([nY, 1]);

% Distribute the calculated values to respective objects with appropriate weighting
for i = 1:N_obj
    y0(id_obj{i}) = y0(id_obj{i}) + w_obj(i) / n_c(c_obj(i)) * y{i};
end

% Calculate the product of matrix A' and vector y
x = single(Dij' * y0);
end
