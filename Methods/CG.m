function [x, k] = CG(A, b, x0, tol, max_iter)
    % A - Coefficient matrix
    % b - Constant vector
    % x0 - Initial guess
    % tol - Tolerance (stopping criterion)
    % max_iter - Maximum number of iterations

    % Initialization
    x = x0;  % Initial solution
    r = b - A*x;  % Residual
    p = r;  % Initial search direction
    rsold = r' * r;  % Initial residual squared norm
    k = 0;  % Iteration count

    for k = 1:max_iter
        Ap = A * p;  % Matrix-vector multiplication A*p
        alpha = rsold / (p' * Ap);  % Compute step size
        x = x + alpha * p;  % Update solution

        r = r - alpha * Ap;  % Update residual
        rsnew = r' * r;  % Compute new residual squared norm

        % Check if the tolerance is satisfied
        if sqrt(rsnew) < tol
            break;
        end

        p = r + (rsnew / rsold) * p;  % Update search direction
        rsold = rsnew;  % Update the residual squared norm
    end
end
