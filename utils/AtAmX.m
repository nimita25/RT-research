function Y = AtAmX(X, AmX, AtmX, var)
% AtAmX Calculates the result of A' * A * X.
%
% INPUT:
%   X     - Input vector.
%   AmX   - Function handle to calculate A * X.
%   AtmX  - Function handle to calculate A' * X.
%   var   - Structure containing relevant information.
%
% OUTPUT:
%   Y - Result of A' * A * X.

% Calculate A * X using the AmX function
AX = AmX(X, var);

% Calculate A' * (A * X) using the AtmX function
Y = AtmX(AX, var);

% Note: Uncomment the following lines if specific conditions apply
% if var.isC == 1 && isreal(Y)
%     Y = complex(Y, zeros(size(Y), 'single'));
% end
end
