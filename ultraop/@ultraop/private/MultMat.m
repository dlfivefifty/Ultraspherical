function M = MultMat(a,bn,varargin)
% MULTMAT(A,BN) This forms a truncation of the Multiplication operator for
% a(x) in Cheb T series.  It truncates so that it works for vectors of length BN. It passes
% back a rectangular matrix. 
% MULTMAT(A,BN,lam) This forms a truncation of the Multiplication operator
% for a(x) in C^{lam} series. 

% not formed. 
    % Computes the multiplication matrix, staying in Chebyshev T. Compute
    % coefficients for a(x)*ChebTfunc
    
    a = a(:); 
    an = length(a);
    %remove trailing zeros
    atrun = truncate(a,eps);
    at = length(atrun); 
    
    %prolong coeffs. 
    %a(length(a)+1:(an+bn)) = 0;
    if(an<bn)
        a(length(a)+1:bn)=0;
    end
    a = a(:);
    if(nargin>2), lam = varargin{1}; end
    
    if(nargin==2||lam==0)
    
        if(bn<=6000)   
        %compute multiplication. 
        if(length(a)>1)
            a=.5*a;
            r=[2*a(1);a(2:end)];
            M = mytoeplitz(r,r);
            a = truncate(a,eps);
            H = hankel(a(2:end));
            %M(2:end,1:end-1) = M(2:end,1:end-1) + H;
            M(2:length(a),1:length(a)-1) = M(2:length(a),1:length(a)-1)+ H;
        else
            M = toeplitz(2*a);
        end
        %half because of J map. 
        %M = M(:,1:bn);
        %MultMat = MultMat(:,1:bn);
        else
    
        %compute multiplication. 
        if(at>1),
            %a=.5*atrun;
            a = .5*a;
            M = sptoeplitz([2*a(1);a(2:end)],[2*a(1);a(2:end)]);
            %Construct hankel by hand for speed. 
            H = sphankel(a(2:end));
            M(2:length(a),1:length(a)-1) = M(2:length(a),1:length(a)-1)+ H;
        else
            M = sptoeplitz(2*atrun);
        end
        %half because of J map. 
        %M = .5*M(:,1:bn);
        end
    elseif(lam==1)
        % Want the U*U Cheb Multiplication matrix.
        % Convert ChebT of a to ChebU 
        a = transMat(length(a),0)*a;
        er = flipud(cumsum(flipud(a(1:2:end))));
        or = flipud(cumsum(flipud(a(2:2:end))));

        r(1:2:length(a)) = er; r(2:2:length(a)) = or; 

        c = [r(3:end) 0 0];
        M = sptoeplitz(r,r) - sphankel(c); 
        M = M(:,1:bn);
    elseif(lam>1)
        % Want the C^{lam}C^{lam} Cheb Multiplication matrix.
        
        % Convert ChebT of a to ChebC^{lam}
        for ind=0:lam-1
            a = transMat(length(a),ind)*a;
        end

        nnza = find(abs(a)>1e-16);
        band = nnza(end);
        M = spalloc(length(a),bn,2*band*bn);
        n=length(a);
        for d=0:n-1
            for k=max(0,d-band):min(n-1,d+band)
                if(k<=d)
                    as=0;
                    c = ChebC(lam,k,d-k,0);
                    for s=0:k
                        if(2*s+d-k<n)
                            as = as + a(2*s+d-k+1)*c;
                            c = c*factorcheb(lam,k,2*s+d-k,s);
                        end
                    end
                    M(d+1,k+1) =  as;
                elseif(k>d)
                   as=0;
                   c = ChebC(lam,k,k-d,k-d);
                   for s=k-d:k
                       if(2*s+d-k<n)
                            as = as + a(2*s+d-k+1)*c;
                            c = c*factorcheb(lam,k,2*s+d-k,s);
                       end
                   end
                   M(d+1,k+1) = as;
                end
            end
        end
        M = M(:,1:bn);
    end 
end