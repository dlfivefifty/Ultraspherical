function fac = factorcheb(v,i,j,s)
%     %c_{s+1}(i,j) = fac*c_{s}(i,j), hopefully can be used to speed up things.
%     fac = ((v+s)/(s+1))*((i-s)/(v+i-s-1))*((2*v+i+j-2*s-2)/(v+i+j-2*s-2));
%     fac = fac*((2*v+i+j-2*s-1)./(v+i+j-2*s-1))*((v+i+j-s-1)/(2*v+i+j-s-1));
%     fac = fac*((j-s)/(v+j-s-1))*((v+i+j-2*s-2)/(i+j-2*s-1))*((v+i+j-2*s-1)/(i+j-2*s));
%     fac = fac*((i+j+v-2*s-2)/(i+j+v-s-1))*((i+j+v-s)/(i+j+v-2*s));
    
      %c_{s+1}(i,j) = fac*c_{s}(i,j), hopefully can be used to speed up things.
      fac = ((v+s)/(s+1))*((i-s)/(v+i-s-1))*((2*v+i+j-s)/(v+i+j-s));
      fac = fac*((v+j-s)/(j-s+1))*((i+j+v-s)/(i+j+v-s+1));
    
end