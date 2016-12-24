function [x nIter]  = newtonMin (f, g, hess, x0, tol, nMax, search)
%
% NEWTONMIN Computes the minimum of a function f by newtons method.
%   [x nIter]  = newtonMin (g, hess, x0, tol, nMax) computes the minimum
%   of f and store each iterations in x and the number of iterations to
%   nIter.
% 

x=[x0];
tolCondition = true;
nIter = 0;
% uses tolCondition so that the loop is done at least once even though x 
% has size 1 < 2
while tolCondition && nIter < nMax
    % computes direction to g
    gk= g(x(:,end));
    p=hess(x(:,end))\-gk;
    
    %computes step size 
    if nargin > 6 && search == true
        % if search=='search'then do weak line search
        a  = weakLineSearch (f, gk, p, x(:,end), 1e-4);
    else
        a = 1;
    end
    
    % update x 
    x(:,end+1)=x(:,end) + a*p;
    
    %update the stop conditions
    tolCondition = norm(x(:,end)-x(:,end-1)) > tol * (1 + norm(x(:,end-1)));
    nIter = nIter + 1;
end
end

