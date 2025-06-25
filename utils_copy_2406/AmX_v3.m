function y = AmX_v3(x, var)
% AmX_v3 Calculates the result of A * X.
%
% INPUT:
%   x   - Beam intensity vector.
%   var - Structure containing information related to A.
%
% OUTPUT:
%   y - Result of the product of A and x.

% Extract necessary information from the input structure
Dij = var.Dij{1};
id_obj = var.id_obj;
N_obj = var.N_obj;

% Calculate the product of matrix A and vector x
y0 = Dij * double(x);

% Initialize a cell array to store results for each object
y = cell(N_obj, 1);

% Distribute the calculated values to respective objects
for i = 1:N_obj
    y{i} = y0(id_obj{i});
end
end
