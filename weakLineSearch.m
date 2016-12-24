function [a]  = weakLineSearch (f, gk, p, x, sigma)
%
% WEAKLINESEARCH finds a step size alpha decent wnough (weak line search)
%   [a]  = weakLineSearch (f, p, sigma) finds a decent stepsize a knowing
%   the function f, the value f the gradient at the current iteration gk,
%   the direction p, the current point x, and a constant sigma used for 
%   Wolfe condition.
%

phi = @(a) f(x + a*p);
% phi'(0)
phiP0=p'*gk;
a = 1;
while phi(a) > phi(0) + sigma*a*phiP0 
    % computes mu such that a(k+1) = mu * a
	mu = -phiP0*a/(2*(phi(a)-phi(0)-a*phiP0));
    % if mu is too small then use backtracking
    if mu < 0.1
        mu = 0.5;
    end
    a = a * mu;
    
    % not good to take to small steps (because will do many iterations for
    % "nothing". Not here it should be better to uses Wolfe's condition
    % for the minimum alpha too, but it's relatively expensive to compute
    % and not so important (in the sense that a very small alpha is still
    % bad)
    if a < 1e-4
        error('WeakLineSearch: Step size too small');
        break;
    end
        
end
 
end

