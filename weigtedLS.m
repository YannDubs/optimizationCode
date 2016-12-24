function [x]  = weigtedLS (A, W, b, method)
%
% WEIGTEDHLS solves a weighted least squares problem .
%   [x]  = weigtedLS (A, W, b, method) computes the weigted least
%   squares solution with QR cholesky 'QRchol 'or QR eigenvalue decomposition 'QReig'.
%   [x]  = weigtedLS (A, W, b) computes the weigted least
%   squares solution with normal equations.
% 

% solves via normal equations
if nargin == 3 
    x = (A'*W*A)\(A'*W*b);
elseif nargin == 4 && strcmp(method,'QRchol') 
	% solves the problem via QR factorization of the reformulated problem.
    % B found by cholesky
    B = chol(W);
    d = B*b;
    C = B*A;
    x = (C'*C)\(C'*d);
elseif nargin == 4 && strcmp(method,'QReig') 
	% solves the problem via QR factorization of the reformulated problem.
    % B found by eigenvalue decomposition
    [U,D] = eigs(W);
    B = U*(D.^(1/2))*U';
    d = B*b;
    C = B*A;
    x = (C'*C)\(C'*d);
else
     error('Missuse of argument: see help weigtedLS');   
end

end

