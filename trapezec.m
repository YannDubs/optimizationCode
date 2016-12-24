function [ intmp ] = trapezec( a, b, M, f )
%MIDPOINT Approximates the integral with formule composite du trapeze
%   Approximates the integral of f between a and b
%   with formule composite du trapeze using M intervalls
H = (b-a)/(M);
Hn = [a:H:b];
sum = 0;
for i = 1:M
    y = (f(Hn(i))+f(Hn(i+1)))/2;
    sum = sum + y;
end
intmp = H*sum;

end
