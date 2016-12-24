function [x maxA k] = cg(x,tol,A,b)
r= b-A*x;
d=r' * r;
bd=b'*b;
k=0;
p=r;
maxA = -inf;

while d > tol^2 * bd
    dOld=d;
    s =A*p; %only multiplication
    a = d./(p'*s);
    x = x + a*p;
    r = r -a*s;
    d = r'*r;
    p = r + (d/dOld) *p;
    k = k+1;
    
    if a > maxA
        maxA = a;
    end
end
end