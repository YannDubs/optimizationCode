function [x, r, n, inc] = newton( F, dF, x0, tol, nmax )
% NEWTON Find the zeros of a of non?linear equation.
% [X] = NEWTON(F,DF,X0,TOL,NMAX) tries to find the zero X of the
% continuous and differentiable function F nearest to X0 using
% the Newton method. DF is a function which take X and return the derivative of F.
% If the search fails an error message is displayed.
%
% [X,R,N,INC] = NEWTON(F,DF,X0,TOL,NMAX) returns the value of the
% residual R in X,the number of iterations N required for computing X and
% INC the increments computed by Newton.
n = 0;
x(n+1) = x0;
r = F(x0);
diff = tol + 1;
while ( diff > tol && n < nmax )
inc(n+1) = F(x(n+1)) / dF(x(n+1));
x(n+2) = x(n+1) - inc(n+1);
diff = abs(inc(n+1));
r = abs(F(x(n+2)));
n = n + 1;
end
if n >= nmax
fprintf(['Newton stopped without converging to the desired tolerance '...
'because the maximum number of iterations was reached\n']);
end
return