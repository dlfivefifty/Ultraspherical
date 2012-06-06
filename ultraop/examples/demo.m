function demo
n=10; %length(a);
lam = 3; 
A = zeros(1); A=sparse(A);
%a = fliplr(chebpoly(chebfun(@(x) sin(100*x)))); %[0 0 1 0 0 0 0 0 0 0]';
%n=length(a); a=a(:);
a = [0 1 zeros(1,98)]';
a=rand(n,1);

nnza = find(abs(a)>1e-16);
band = nnza(end);

for d=0:n-1
    for k=max(0,d-band):min(n-1,d+band)
        if(k<=d)
            as=0;
            c = ChebC(lam,k,d-k,0);
            for s=0:k
                if(2*s+d-k<n)
                    as = as + a(2*s+d-k+1)*c;%ChebC(lam,k,2*s+d-k,s);
                    c = c*factorcheb(lam,k,2*s+d-k,s);
                end
            end
            A(d+1,k+1) =  as;
        elseif(k>d)
           as=0;
           c = ChebC(lam,k,k-d,k-d);
           for s=k-d:k
               if(2*s+d-k<n)
                    as = as + a(2*s+d-k+1)*c;%ChebC(lam,k,2*s+d-k,s);
                    c = c*factorcheb(lam,k,2*s+d-k,s);
               end
           end
           A(d+1,k+1) = as;
        end
    end
end
B = transMat(2*n,2)*transMat(2*n,1)*transMat(2*n,0)*full(MultMat(transMat(n,0)\(transMat(n,1)\(transMat(n,2)\a)),2*n,0))*inv(transMat(2*n,0))*inv(transMat(2*n,1))*inv(transMat(2*n,2));
%B = transMat(100,0)*full(MultMat(transMat(10,0)\a,100,0))*inv(transMat(100,0));
%A=A*diag(.25*ones(n,1)); 
%spy(isnan(A))
norm(A-B(1:n,1:n))
end



