function spy(N,n)
% SPY plot
%
% spy(N), spy plot of a 100 by 100 discretisation of an ultraop. 
%
% spy(N,n), spy plot of a n by n discretisation of an ultraop. 

ish = ishold; 
if nargin < 2, n = 100; end % default to n=100 realisation. 

spy(realisation(N,n)); 

if ~ish, hold off; end 
% if nargout > 0, varargout = {h}; end
end