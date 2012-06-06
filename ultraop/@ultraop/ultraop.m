% ULTRAOP Constructor for ultraop
% 
% ULTRAOP(F) constructs a ultraop object. 

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

