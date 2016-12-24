function [ intmp ] = midpointc( a, b, M, f )
%MIDPOINT Approximates the integral with formule composite du point milieu
%   Approximates the integral of f between a and b
%   with formule composite du point milieu using M intervalls
H = (b-a)/(M);
Hn = [a:H:b];
sum = 0;
for i = 1:M
    x = (Hn(i)+Hn(i+1))/2;
    
    sum = sum + f(x);
end
intmp = H*sum;

end


