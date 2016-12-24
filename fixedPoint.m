function [p, res, niter] = fixedPoint ( phi, x0, tol, nmax )
% FIXEDPOINT Fixed point iterations.
% P = FIXEDPOINT( PHI, X0, TOL, NMAX ) tries to find the fixed point P of the
% continuous and differentiable function PHI using fixed point iterations
% starting from X0. PHI accepts real scalar input x and returns a real
% scalar value. If the search fails an errore message is displayed.
%
% [P, RES, NITER] = FIXEDPOINT(PHI,...) returns the value of the residual in P
% and the iteration number at which P was computed.
%
p = x0;
niter = 0;
res = tol + 1;
while (res > tol )
    p = phi (p) ;
    res = abs(phi (p) - p);
    niter = niter + 1;
    
    if niter >= nmax
        fprintf ('The fixedPoint didn''t converge in %1.1f iterations \n', nmax);
        break;
    end
    
end
return;