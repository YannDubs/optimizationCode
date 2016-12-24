function [t, y] = heun( fun, tspan, y0, N )
%  HEUN Solve differential equations using the HEUN method.
%  [T, Y] = HEUN( FUN, TSPAN, Y0, N ), with TSPAN = [T0, TF],
%  integrates the system of differential equations y'=f(t, y) from time T0 
%  to time TF, with initial condition Y0, using the HEUN method on 
%  an equispaced grid of N intervals. Function FUN(T, Y) must return
%  a column vector corresponding to f(t, y). Each row in the solution 
%  array Y corresponds to a time returned in the column vector T.
%  [T, Y] = HEUN( FUN, TSPAN, Y0, N, P1, P2,...) passes the additional
%  parameters P1, P2, ... to the function FUN as FUN(T, Y, P1, P2, ...).
%
%  NOTE: 
%  if for example
%      fun = @(t,y) sin(t*y) + y;
%  then
%      fun(a,b);
%  evaluates the function fun at the points t = a, y = b.

% time step
h = ( tspan(2) - tspan(1) ) / N;

% time snapshots 
t = linspace( tspan(1), tspan(2), N+1 );

% initialize the solution vector
y = [y0 zeros(1,N)];

% time loop (n=0,...,n, but array indeces in Matlab start at 1)
for n = 1:N
    y(n+1) = y(n) + 0.5 * h * ( fun( t(n), y(n) ) + fun( t(n)+h, y(n) + h*fun( t(n), y(n) ) ) );
end

return