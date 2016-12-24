function [ t u ] = backwardEulerSystem( F,tspan,y0,N )
%BACKWARDSYSTEM Summary of this function goes here
%   Detailed explanation goes here
h = abs(tspan(2)-tspan(1))/N;
t = linspace(tspan(1),tspan(2),N+1);
u=[y0];
options.Display='off';
options.TolFun = 1e-9;
options.MaxFunEvals = 1e4;
for n=1:N
    f = @(x) -x+h.*F(t(n+1),x)+u(:,n);
    u(:,n+1)= fsolve (f,u(:,n),options);
end

end


