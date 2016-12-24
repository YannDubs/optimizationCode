function [x nIter] = bfgs(f, g, x0, tol, nMax, search)
%
% BFGS Computes the minimum of a function f by bfgs method.
%   bfgs(f, x0, tol, nMax,search) computes the minimum of f and
%   stores each iterations in x and the number of iterations to
%   nIter.
%

x = [x0];
tolCondition = true;
nIter = 0;
% gx contains the value of g(x(k+1) and g(x(k))
% simply for initializing: x0 will never be used as gx 
gx = [g(x0) x0];
[~, n] = size(x0);
C = eye(n); 

while tolCondition && nIter < nMax 
    % computes direction to g
    p = -C*gx(:,1);
    
    %computes step size 
    if nargin > 5 && search == true
        % if search=='search'then do weak line search
        a  = weakLineSearch (f, gx(:,1), p, x(:,end), 1e-4);
    else
        a = 1;
    end
    
    % update x 
    s = a*p;
    x(:,end+1) = x(:,end) + s;
    
    % updates the gradient values gx
    gx(:,2) = gx(:,1);
    gx(:,1) = g(x(:,end));
    
    % update the inverse of the matrix
    y = gx(:,1)-gx(:,2);
    rho = 1/(y'*s);
    C = (eye(n)-rho*s*y')*C*(eye(n)-rho*y*s') + rho*s*s';
    
    % update the stop conditions
    tolCondition = norm(x(:,end)-x(:,end-1)) > tol * (1 + norm(x(:,end-1)));
    nIter = nIter + 1;
end
end