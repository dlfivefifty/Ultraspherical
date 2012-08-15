function [Wa Wb r] = recoverCoeffsBC(L)
%RECOVERCOEFFSBC  Recover coefficient functions of a linear BCs
% [WA WB R] = RECOVERCOEFFS(L) returns, for a linear operator L on a domain 
% [a,b], two matrices WA and WB such that
%   sum_{k=0}^{n} ( WA(:,k)*u^(k)(a) + WB(:,k)*u^(k)(b) ) = R.
%
% Example 1:
%   N = chebop(@(x,u) diff(u,2),[-1 1]);
%   N.lbc = @(u) u-2;
%   N.rbc = @(u) 3*diff(u)-4;
%   [Wa Wb r] = recoverCoeffsBC(N)
%   Wa =                Wb =                r = 
%          1     0             0     0            2
%          0     0             0     3            4
%
% Example 2:
%   N = chebop(@(x,u) diff(u,2),[-1 1]);
%   N.lbc = @(u) diff(u) - 1;
%   N.bc = @(u) u(1) - u(-1);
%   [Wa Wb r] = recoverCoeffsBC(N)
%   Wa =                Wb =                r = 
%          0     1             0     0            1
%         -1     0             1     0            0

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

if isa(L,'chebop'), L = linop(L); end
if isa(L,'ultraop')
   L1 = chebop(@(x,u) diff(u,L.order)); 
   L1.lbc = L.lbc; L1.rbc=L.rbc;
   L = linop(L1); 
end


dom = L.fundomain;
m = L.difforder;

N = 3*m;
[LL BC] = feval(L,N);

B = zeros(2*m,size(LL,2));
for k = 0:m-1
    Dk = feval(diff(dom,k),N);
    B(2*k+1,:) = Dk(1,:);
    B(2*k+2,:) = Dk(end,:);
end

numbc = size(BC,1);
valid = false(numbc,1);
for k = 1:numbc
    valid(k) = rank([B ; BC(k,:)]) ~= 2*m+1;
end

N = 2*m;
[LL BC c] = feval(L,N);
B = zeros(2*m,size(LL,2));
for k = 0:m-1
    Dk = feval(diff(dom,k),N);
    B(2*k+1,:) = Dk(1,:);
    B(2*k+2,:) = Dk(end,:);
end

coeffs = BC(valid,:)/B;
coeffs(abs(coeffs)<1e-14) = 0;

r = c(valid);
Wa = coeffs(:,1:2:end);
Wb = coeffs(:,2:2:end);

end