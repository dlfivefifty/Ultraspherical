function spy(N,n)
% spy plot of a ultraop. This is the spy plot of a n by n realisation.

ish = ishold; 
if nargin < 2, n = 100; end % default to n=100 realisation. 

spy(realisation(N,n)); 

if ~ish, hold off; end 
% if nargout > 0, varargout = {h}; end
end