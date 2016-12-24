function [x,inc_vec,iter] = newtonSys(F, J, x0, tol, maxiter)
% NEWTONSYS Find the zeros of a nonlinear equations.
% [X] = NEWTONSYS(F,J,X0,TOL,NMAX) tries to find the zero X of the
% system of continuous and differentiable functions F nearest to X0 using
% the Newton method. J is a function which takes X and return the Jacobian
% of F evaluated at X. If the method fails an error message is displayed.
%
% X0 is a column vector of dimension n x 1
%
% F is a function handle F = @(x)[F_1(x);
% F_2(x);
% ...;
% F_n(x)];
% such that F(x) is a vector of size n x 1
%
% J is a function handle J = @(x)[J_11(x) ... J_1n(x);
% ... ... ... ;
% J_n1(x) ... J_nn(x)];
% such that J(x) is a matrix of size n x n
%
% [X,INC_VEC,ITER] = NEWTONSYS(F,DF,X0,TOL,NMAX) returns a vector INC_VEC
% that contains the norm of the increments and the the number of
% iterations ITER required for computing X.
inc_vec = [];
iter = 0;
x_old = x0;
x = x0;
norm_increment = tol + 1;
while ( norm_increment > tol && iter < maxiter )
iter = iter + 1;
r = F(x_old);
increment = -J(x_old)\r;
x = x_old + increment;
norm_increment = norm(increment);
inc_vec = [inc_vec, norm_increment];
x_old = x;
end
if iter > maxiter
fprintf(['Newton stopped without converging to the desired tolerance '...
'because the maximum number of iterations was reached\n']);
end
end