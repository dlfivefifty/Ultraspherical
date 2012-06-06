function diffMat = diffMat(n,lam)
%DIFFMAT(N,LAM) This computes the truncation of the differentiation operator
% between C^{lam} to C^{lam+1}. It returns an n x n matrix. 

    if(n==0)
        diffMat=[];
        return;
    end
    if(n<1e2)
        if(lam==0)
            diffMat = diag((1:n-1)',1);
        else
            diffMat = diag(2*lam*ones(n-1,1),1); 
        end
    else
    %Cheb T is a special case.
        if(lam==0)
            diffMat = spdiags((0:n-1)',1,n,n); 
        else
            diffMat = spdiags(2*lam*ones(n,1),1,n,n); 
        end
    
    end
end