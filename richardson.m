function [x, relres, iter] = richardson ( A, b, x0, alpha1, tol, maxiter )
% output: x, solution
% relres, vector of relative residuals norm(r_k)/norm(r_0)
% iter, number of iterations
%
% intput: A, matrix
% b, right hand side
% x0, initial guess
% alpha, acceleration parameter
% tol, tolerance
r0 = ( b - A*x0 ); % initial residual
relres = [];
relres = [relres, norm(r0)/norm(r0)];
iter = 1;
x = x0;
residual = r0;
while (relres(end) > tol && iter < maxiter)
    x = x + alpha1 * residual; 
    residual = (b - A*x); % residual
    relres = [relres, norm(residual)/norm(r0)]; % relative residual
    iter = iter + 1;
end
end
