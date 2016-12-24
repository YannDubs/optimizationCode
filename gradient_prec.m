function [x,relres,iter] = gradient_prec ( A, b, x0, tol, maxiter, P )
% output: x, solution
% relres, vector of relative residuals norm(r_k)/(r_0)
% iter, number of iterations
%
% intput: A, matrix
% b, right hand side
% x0, initial guess
% tol, tolerance
% P, preconditioner
r0 = b - A*x0;
relres = [];
relres = [relres, norm(r0)/norm(r0)];
iter = 1;
x = zeros(size(x0));
r = r0;
while (relres(end) >tol && iter < maxiter)
z = P\r;
alpha = (z'*r)/(z'*(A*z));
x = x + alpha*z;
r = r -alpha*A*z;
relres = [relres, norm(b-A*x)/norm(r0)]; % relative residual
iter = iter + 1;
end
end