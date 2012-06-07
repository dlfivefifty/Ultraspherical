function c3 = ChebC(v,m,n,s)
        
        %algebraic mess:
        c3 = ((m+n+v-2*s)/(m+n+v-s))*prod(((v:v+s-1)./(1:s)).*((2*v+m+n-2*s:2*v+m+n-s-1)./(v+m+n-2*s:v+m+n-s-1)));
        c3 = c3*prod(((v:v+m-s-1)./(1:m-s)).*((n-s+1:m+n-2*s)./(v+n-s:v+m+n-2*s-1)));       
end