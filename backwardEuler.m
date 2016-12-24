
%{ 
function [t, u] = backwardEuler(fun, tspan, y0, Nh)

%BACKWARDEULER Solve differential equations using the backward
% Euler method and non?linear solves by the fixpoint iterations
%
% [T,U]=BACKWARDEULER(FUN,TSPAN,Y0,NH) with TSPAN=[T0,TF]
% integrates the system of differential equations
% y'=f(t,y) from time T0 to TF with initial condition
% Y0 using the backward Euler method on an equispaced
% grid of NH intervals.Function FUN(T,Y) must return
% a column vector corresponding to f(t,y). Each row in
% the solution array Y corresponds to a time returned
% in the column vector T.
%
% note:
% if for example
% fun = @(t,y) sin(t*y) + y;
% then
% fun(a,b)
% evaluate the function odefun at the points t=a, y=b.
% time step
% note: n=0,...,Nh, but array indeces in Matlab start at 1!
h=( tspan(2)-tspan(1) ) / Nh;
% time snapshots
t=linspace(tspan(1),tspan(2),Nh+1);
% define solution vector
u(1)=y0;
% parameters
tol = 1e-9;
nmax = 10000;
% time loop
for n = 1:Nh
phi = @(x) ( u(n) + h * fun(t(n+1),x) );
u(n+1) = fixedPoint(phi,u(n),tol,nmax);
end
% finish
return
%}


function [t, y] = backwardEuler( fun, tspan, y0, N )
% BACKWARDEULER Solve differential equations using the backward Euler method.
% [T, Y] = BACKWARDEULER( FUN, TSPAN, Y0, N ), with TSPAN = [T0, TF],
% integrates the system of differential equations y'=f(t, y) from time T0
% to time TF, with initial condition Y0, using the backward Euler method on
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

h=(tspan(2)-tspan(1))/N;
t=linspace(tspan(1), tspan(2), N+1);
u(1)=y0;

options=optimset;
options.Display='off';
options.TolFun=1.e-09;
options.MaxFunEvals=10000;
    
for n=1:N
    F = @(x) x-u(n)-h*fun(t(n+1),x);
    u(n+1) = fsolve(F,u(n) + h * fun( t(n), u(n) ),options);
end
y=u;





