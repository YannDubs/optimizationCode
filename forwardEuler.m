function [t, y] = forwardEuler( fun, tspan, y0, N )
% FORWARDEULER Solve differential equations using the forward Euler method.
% [T, Y] = FORWARDEULER( FUN, TSPAN, Y0, N ), with TSPAN = [T0, TF],
% integrates the system of differential equations y'=f(t, y) from time T0
% to time TF, with initial condition Y0, using the forward Euler method on
% an equispaced grid of N intervals. Function FUN(T, Y) must return
% a column vector corresponding to f(t, y). Each row in the solution
% array Y corresponds to a time returned in the column vector T.
%
% NOTE:
% if for example
% fun = @(t,y) sin(t*y) + y;
% then
% fun(a,b);
% evaluates the function fun at the points t = a, y = b.

t=linspace(tspan(1),tspan(2),N+1);
h=(tspan(2)-tspan(1))/N;
u(1)=y0;
i=1;
for n = 1:N
    u(n+1) = u(n) + h * fun( t(n), u(n) );
end
y=u;
