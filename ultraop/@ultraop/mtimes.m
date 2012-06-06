function f = mtimes(N,f)
% * mtimes, forward action of an ultraop
%
% N*f returns the chebfun resulting from applying the operator N to the
% chebfun f. 

if(~isa(f,'chebfun'))
    error('Ultraop:mtimes','Not a chebfun\n');
else
    cf=fliplr(chebpoly(yf));
    lf=length(cf);
    %cy = prolong(cy,ly);
    f = realisation(N,lf)*cf';
    %y = truncate(y,eps);
    for s = N.order-1:-1:0
        f=transMat(lf,s)\f;
    end
    f = chebfun(flipud(f),'coeffs');
end
end