function [ t,u ] = ForwardEulerSystem( F, tspan, y0, N )
%FORWARDEULERSYSTEM Summary of this function goes here
%   Detailed explanation goes here

t = linspace (tspan(1), tspan(2), N+1);
h = (tspan(2) - tspan(1))/N;
u(:,1)=y0;
for n = 1:N
    u(:,n+1) = u(:,n) + h.*F(n*h,u(n));
end

end

