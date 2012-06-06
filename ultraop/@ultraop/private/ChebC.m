function c3 = ChebC(v,m,n,s)
        if(m+n-2*s<0 || m-s<0 || n-s<0), c=0; return; end
%         %oldest way 
        %c = (m+n+v-2*s)*Poch(v,s)*Poch(v,m-s)*Poch(v,n-s)*Poch(2*v,m+n-s)*factorial(m+n-2*s);
    	%c3= c/((m+n+v-s)*factorial(s)*factorial(m-s)*factorial(n-s)*Poch(v,m+n-s)*Poch(2*v,m+n-2*s));
%         
%         %older way
         %c1 = (m+n+v-2*s)*Poch(v,s)*Poch(v,m-s)*prod(2*v+m+n-2*s:2*v+m+n-s-1)*factorial(m+n-2*s);
     	%c1= c1/((m+n+v-s)*factorial(s)*factorial(m-s)*factorial(n-s)*prod(v+n-s:v+m+n-s-1))
%         
%         %better way
        %c2 = (m+n+v-2*s)*Poch(v,s)*Poch(v,m-s)*prod(2*v+m+n-2*s:2*v+m+n-s-1)*prod(n-s+1:m+n-2*s);
     	%c2= c2/((m+n+v-s)*factorial(s)*factorial(m-s)*prod(v+n-s:v+m+n-s-1))
%         
%         %even better way.
         %c3 = ((m+n+v-2*s)/(m+n+v-s))*(Poch(v,s)/factorial(s));
         %c3 = c3*(Poch(v,m-s)/factorial(m-s))*((prod(2*v+m+n-2*s:2*v+m+n-s-1))/prod(v+n-s:v+m+n-s-1))*prod(n-s+1:m+n-2*s);
         %abs(c2-c3);
%         
%         %no overflow way. 
         %c3 = ((m+n+v-2*s)/(m+n+v-s));
         %c3 = c3*prod((v:v+s-1)./(1:s));
         %c3 = c3*prod((v:v+m-s-1)./(1:m-s)); 
         %c3 = c3*prod((2*v+m+n-2*s:2*v+m+n-s-1)./(v+m+n-2*s:v+m+n-s-1));
         %c3 = c3*prod((n-s+1:m+n-2*s)./(v+n-s:v+m+n-2*s-1))
        
        %no overflow and faster
        c3 = ((m+n+v-2*s)/(m+n+v-s))*prod(((v:v+s-1)./(1:s)).*((2*v+m+n-2*s:2*v+m+n-s-1)./(v+m+n-2*s:v+m+n-s-1)));
        c3 = c3*prod(((v:v+m-s-1)./(1:m-s)).*((n-s+1:m+n-2*s)./(v+n-s:v+m+n-2*s-1)));
        %c3 = nchoosek(m+n-2*s,n-s)./factorial(s)*c3;
        
end

function p = Poch(v,k)

 p = prod(v:v+k-1);
end