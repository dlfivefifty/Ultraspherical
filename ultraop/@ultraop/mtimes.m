function y = mtimes(N,y)
% Computes N*y.
if(~isa(y,'chebfun'))
    error('Ultraop:mtimes','Not a chebfun\n');
else
    cy=fliplr(chebpoly(y));
    ly=length(cy);
    %cy = prolong(cy,ly);
    y = realisation(N,ly)*cy';
    %y = truncate(y,eps);
    for s = N.order-1:-1:0
        y=transMat(ly,s)\y;
    end
    y = chebfun(flipud(y),'coeffs');
end
end