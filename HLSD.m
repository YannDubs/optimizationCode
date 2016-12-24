function [x maxA k] = HLSD(x,tol,A,b)
r= b-A*x;
d=r' * r;
bd=b'*b;
k=0;
maxA = -inf;

while d > tol^2 * bd
    s =A*r;
    
    if mod(k,2) == 0 %if even update if not keep
        a = d./(r'*s); 
    end
    
    x = x + a*r;
    r = r -a*s;
    d = r'*r;
    k = k+1;
    
    if a > maxA
        maxA = a;
    end
end
end