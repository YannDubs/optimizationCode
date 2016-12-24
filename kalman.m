function [x P]  = kalman (x0, P0, Gk, Hk, z, Qk, Rk, sk, method)
[ nVar, ~ ]=size(Gk);
[ ~, N ] = size(z);
x = [x0];
P = [P0];
I=eye(nVar);
%% Original Kalman Algorithm
if nargin == 8
    for k=1:N
        % computes predicted estimates
        xkp = Gk * x(:,k) + sk ;
        % computes predicted covaraince
        Pkp = Gk * P(:,:,k) * Gk' + Qk ;
        wk = z(k) - Hk * xkp;
        Sk = Hk * Pkp * Hk' + Rk;
        Kk = Pkp*Hk'*inv(Sk);
        x(:,k+1) = xkp + Kk*wk;
        P(:,:,k+1) = (I-Kk*Hk)*Pkp;
    end     
%% Solves weighted Least Squares Problem
%Note: the function could have been writtern in a nicer way(in the same way 
% as for calling block_trid) but it would have put in memory every A and 
% B but this is unnecessary. Here I do not store any unnecessary array and do a single loop.
elseif nargin == 9 && strcmp(method,'WLS') 
    % Note: k_1 represents k+1
    for k=1:N        
        % Note that Rk and Qk and Gk doesn't change at each step here so constant
        % value
        % we compute Ak and ck as if it were the last iteration (AN and CN 
        % at each time in order to be able to graph xk)
        Ak = inv(Qk) + Hk'*inv(Rk)*Hk ;
        Bk_1 = -inv(Qk)*Gk;
        ck = Hk'*inv(Rk)*z(k) + inv(Qk)*sk ;
        
        if k==1
            % C1
            ck = ck + inv(Qk)*(Gk*x0);
            %initializes the c and A update as first c and k
            A_update=Ak;
            c_update=ck;
        end

        % This loop goes from k=2 to N (from i = 1: N-1)
        if k~=1
            % here 'in reality k+1 is k', in the sense that nothing was yet
            % updated with k (next loop) so we are working 'with the previous k.'
            Lk = Bk*inv(A_update);
            A_update = Ak - Lk*Bk';
            c_update = ck - Lk*c_update;
        end
        Bk = Bk_1;
        
        P(:,:,k+1)=inv(A_update);
        x(:,k+1)=P(:,:,k+1)*c_update;
        
        % These lines add the term which is different between any k and k=N
        % at the last iteration this is not necessary
        if k~=N
            A_update = A_update + Gk'*inv(Qk)*Gk;
            c_update = c_update - Gk'*inv(Qk)*sk;
        end 
    end

%% Kalman smoother
% Note: uses backward ssubstitution to smooth the function
elseif nargin == 9 && strcmp(method,'KS') 
    for k=1:N-1
        A(:,:,k) = inv(Qk) + Hk'*inv(Rk)*Hk + Gk'*inv(Qk)*Gk;
        B(:,:,k+1)= -inv(Qk)*Gk;
        Bt(:,:,k+1)= (-inv(Qk)*Gk)';
        c(:,k)=Hk'*inv(Rk)*z(k)+inv(Qk)*sk-Gk'*inv(Qk)*sk;
    end
        c(:,1) = c(:,1)+inv(Qk)*(Gk*x0);
        c(:,N) = Hk'*inv(Rk)*z(N)+inv(Qk)*sk;
        A(:,:,N) = inv(Qk) + Hk'*inv(Rk)*Hk;   
        % The function from the first assignement
        x = block_trid(A(:,:,:), B(:,:,2:end), Bt(:,:,2:end), c);
        x = [x0 x];
else
     error('Missuse of argument');   
end

end

