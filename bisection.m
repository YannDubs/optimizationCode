function [zero,res,nit] = bisection(fun,a,b,tol,nmax)
% BISECTION Find a zero of a nonlinear scalar function inside an interval.
%   ZERO=BISECTION(FUN,A,B,TOL,NMAX) finds a zero of the continuous
%   function FUN in the interval [A,B] using the bisection method and returns
%   a number ZERO that is the zero of the function.
%   FUN accepts real scalar input x and returns a real scalar value;
%   FUN can also be an inline object.
%   TOL is the tolerance on error allowed and NMAX the maximum number of iterations.
%   If the search fails an error message is displayed.
%
%   [ZERO,RES,NIT]=BISECTION(FUN,...) returns also RES, that is the residual
%   evaluated at the last iterate, and NIT the number of iterations.
%
%

% input check
if a>=b
    error('WARNING: b doit etre plus grand que a');
end
fa = fun(a);
fb = fun(b);
if (fa * fb > 0)
    error('WARNING: f(a),f(b) doivent avoir signs differents!');
end
nit = 0;
if fa == 0
    zero = a;
    res = 0;
    return
elseif fb == 0
    zero = b;
    res = 0;
    return
end
%initialize error
I = tol +1;
while (I > tol && nit < nmax)
    
    x = 0.5 * (a+b);
    fx = fun(x);
    nit = nit +1;
    zero(nit) = x;
    
    if (fx == 0)
        res = 0;
        return;
    end
    
    if fa * fx < 0
        b = x;
    else
        a = x;
    end
    I = b-a;
end

res = fx;
if nit >= nmax
    error('WARNING: depass? numero maximal d''iterations, covergence n''est pas garantie ');
end

        
end