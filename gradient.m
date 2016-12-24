function [x,relres,iter] = gradient ( A, b, x0, tol, maxiter )
% output: x, solution
% relres, vector of relative residuals norm(residual_k)/(residual_0)
% iter, number of iterations
%
% intput: A, matrix
% b, right hand side
% x0, initial guess
% tol, tolerance
% maxiter, maximum number of iterations

r0 = b - A*x0; % intial residual
relres = [];
relres = [relres, norm(r0)/norm(r0)];
iter = 1;
x = zeros(size(x0));
r = r0;

while (relres(end) >tol && iter < maxiter)
    alpha = (r'*r)./(r'*A*r); % acceleration parameter
    x = x + alpha.*r; % iteration gradient
    r = (b-A*x); % residual
    relres = [relres, norm(r)/norm(r0)]; % relative residual
    iter = iter + 1;
end

end