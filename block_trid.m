function x = block_trid (A,C,B,b)
%
% Computes the LU block factorizationand then the soltion
% of Mx=b. With M a block tridiagonal matrix of m*n and A,C,B  
% the set of m*m matrices which corresponds respectively to
% the center, lower and upper diagonal.
%
[m,m,n] = size(A);
U=[];L=[];

% Computes LU factorization by block | M=LU
U(:,:,1)=[A(:,:,1)];
y=[b(1:m);0];
for i=1:n
    % uses the ifs to use only 2 loops (not 3)
    if i~= n 
        L(:,:,i)=C(:,:,i)/U(:,:,i);
        U(:,:,(i+1))=A(:,:,(i+1))-L(:,:,i)*B(:,:,i);
    end
    
    % computes directly y | Ly=b
    if i~=1
        y((i-1)*m+1:i*m)=b((i-1)*m+1:i*m)-L(:,:,i-1)*y((i-2)*m+1:(i-1)*m);
    end
end

% computes x | Ux=b, and returns the value
x=zeros(m*n,1);
x((n-1)*m+1:n*m)=[U(:,:,n)\y((n-1)*m+1:n*m)];
for i=n-1:-1:1
    x((i-1)*m+1:i*m)=U(:,:,(i))\(y((i-1)*m+1:i*m)-B(:,:,i)*x((i*m+1:(i+1)*m)));
end

end