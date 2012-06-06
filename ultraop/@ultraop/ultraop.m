% ULTRAOP Constructor for ultraop
% 
% ULTRAOP(F) constructs a ultraop object. 
%
% N=ULTRAOP(F), where F is a n by 1 cell array of function handles or 
% chebfuns will set-up an ultraop representing the differential operator
%
%  N(u) = F{1}du^{n}/dx + F{2}du^{n-1}/dx + ... + F{end}u.
%
% N=ultraop(F), where F is a function handle will return an operator
% representing that handle. We follow the same syntax as the chebop object. 
%
% N=ultraop(F), where F is a chebfun/linop or a linear chebfun/chebop will 
% convert a differential operator represented as a linop or chebop into an
% ultraop.  This is, perhaps, easier syntax for users who are used to the 
% chebfun syntax. 
%
% Example of solving differential equation:
%
% N = ultraop(@(x,u) diff(u,2) + x.*u); N.lbc = 1; N.rbc = 1; 
% u = N\0; 
%
% See chebfun/chebop.

classdef ultraop
    properties ( GetAccess = 'public' , SetAccess = 'public' )
        order  %Order of DE
        lbc    %Left boundary condition
        rbc    %Right boundary condition
        DEcoeffs %Coefficients of DE (quasi-matrix of chebfuns)
        num    %Solving DE nonadaptively with num points.
    end
    methods
        function N = ultraop ( varargin )
            if( nargin == 0 )
            else
                N = ctor(N , varargin{:} );  % pass to constructor.
            end
        end
    end %end methods.
end

