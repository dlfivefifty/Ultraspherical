function u = prolong(u,n)
% PROLONG(U,N) This either adds zeros onto u to make it of length n or it
% truncates u to form vector of length n. 

        % prolongs or truncates by adding zeros; 
        u = u(:);
        if(length(u)>=n)
            %restrict vector. 
            u = u(1:n);
        else
             u(length(u)+1:n)=0;
        end
end