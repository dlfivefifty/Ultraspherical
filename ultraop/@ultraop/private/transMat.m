function spM = transMat(n,lam)
% TRANSMAT(N,LAM) This computes the truncation of the operator that 
% transforms C^{lam} (Ultraspherical polynomials) to C^{lam+1}.  The 
% truncation gives back a matrix of size n x n. 

if(n==1)
    spM = 1;
    return;
end

%Relation is C_n^(lam) = (lam/(n+lam))(C_n^(lam+1) - C_{n-2}^(lam+1))

if(n<1e2)
    if(lam~=0)
        dg = lam./(lam + (2:n-1))';
        spM = diag([1;lam./(lam+1);dg],0) + diag(-dg,2);
    elseif(lam==0)
        dg = .5*ones(n-2,1);
        spM = diag([1;.5;dg],0) + diag(-dg,2);
    end
else
if(lam ~= 0) 
    dg = lam./(lam + (2:n-1))'; 
    B(:,1) = [1;lam./(lam+1);dg]; B(:,2)=[0;0;-dg];
    spM = spdiags(B,[0,2],n,n);
elseif(lam==0)
    % Cheb T is special case because of different scaling.
    dg = .5*ones(n-2,1);
    B(:,1) = [1;.5;dg]; B(:,2)=[0;0;-dg];
    spM = spdiags(B,[0,2],n,n);
end
end
end 