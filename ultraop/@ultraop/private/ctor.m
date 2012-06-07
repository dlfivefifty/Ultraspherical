function N = ctor(N,varargin)
% This is the ultraop constructor. It sets up a ultraop object
% representing a DE operator.
 
N.order = 0; N.DEcoeffs = [];
N.lbc=[]; N.rbc=[]; N.num=Inf;
k=2;
if nargin>1
    if(isnumeric(varargin{1}))
        N.num=varargin{1}; k=3;
    end
end
while k<nargin
    if(strcmpi('order',varargin{k-1}))
        N.order = varargin{k};
        k=k+2;
    elseif(strcmpi('lbc',varargin{k-1}))
%         N.lbc = varargin{k};
        N.lbc = tidylbcs(varargin{k});
        k=k+2;
    elseif(strcmpi('rbc',varargin{k-1}))
%         N.rbc = varargin{k};
        N.rbc = tidyrbcs(varargin{k});
        k=k+2;
    else
        k=k+1;
    end
end

DE = varargin{:}; 
if(iscell(DE))
    % cell array of function handles
    for j=1:length(DE)
        DEcheb(:,j) = chebfun(DE{j});
        % compute the true polynomial interpolant.
        DEcheb(:,j) = chebfun(DE{j},length(DEcheb(:,j)));
    end
    N.DEcoeffs = DEcheb; DE = DEcheb;
elseif isa(DE,'function_handle')
    % If it is a function handle, let chebop deal with it and then convert
    % to a ultraop. 
    DE = chebop(DE); 
    if(~islinear(DE))
        error('ULTRAOP:CTOR:INPUT','Chebop must be linear.');
    end
    DE = recoverCoeffs(DE); DE=DE{:}; DE = DE(:,end:-1:1);
%     % collect lbc and rbc. 
%     if ~isempty(DE.lbc)
%         N.lbc = tidylbcs(DE.lbc);
%     end
%     if ~isempty(DE.rbc)
%         N.rbc = tidyrbcs(DE.rbc);
%     end
elseif isa(DE,'linop')
    % construct the ultraop from a linop. 
    DE = recoverCoeffs(DE); DE=DE{:}; DE = DE(:,end:-1:1);
    % collect lbc and rbc. 
    if ~isempty(DE.lbc)
        N.lbc = tidylbcs(DE.lbc);
    end
    if ~isempty(DE.rbc)
        N.rbc = tidyrbcs(DE.rbc);
    end
elseif isa(DE,'chebop')
    % construct the ultraop from a linear chebop. 
    if(~islinear(DE))
        error('ULTRAOP:CTOR:INPUT','Chebop must be linear.');
    else
       DE = linop(DE); 
    end
    DE = recoverCoeffs(DE); DE=DE{:}; DE = DE(:,end:-1:1);
    % collect lbc and rbc. 
    if ~isempty(DE.lbc)
        N.lbc = tidylbcs(DE.lbc);
    end
    if ~isempty(DE.rbc)
        N.rbc = tidyrbcs(DE.rbc);
    end
else
    error('ULTRAOP:CTOR:INPUT','Can only construct a ultraop from a cell array of function handles, linop or chebop');
end

% extra error in case things unexpectly didn't parse.
if(~isa(DE(:,end),'chebfun'))
    error('DE should be a cell array or a quasi-matrix of chebfuns');
end

% Construct the ultraop.
if(N.order==0)
    if(~isempty(DE))
        %User wants order to be computed.
        N.order = size(DE,2)-1;
        N.DEcoeffs = DE;
    else
        %give back empty ultraop.
        return;
    end
end
if(N.order>0)
    % User wants that order, and default (=0) DEcoeffs.
%     chebzero = chebfun(@(x) 0*x);
    %chebzero = @(x) 0*x;
    while(size(DE,2)<N.order)
%         DE =[chebzero DE];
        error('ULTRAOP:CTOR:ORDER','Not enough variable coefficients for DE order.\n');
    end
    if(size(DE,2)>N.order+1)
        %truncate the end of DE (higher order coeffs)
%         DE = DE(:,1:N.order+1);
        error('ULTRAOP:CTOR:ORDER','Too many variable coefficients for DE order.\n');
    end
end
N.DEcoeffs = DE;
end


function [lbcfh] = tidylbcs(lbc)
% assigns the boundary conditions. Returns them as function handle.

if isdouble(lbc)
    if numel(lbc)==1
        lbcfh = @(u) u(-1) - lbc; 
    else
       error('LBC:Input','Left boundary condition is not of robin-type'); 
    end
elseif isa(lbc,'function_handle')
    lbcfh = lbc; % assume that it is in the correct format. 
else
   error('LBC:Input','Left boundary condition is not function handle or double'); 
end

end

function [rbcfh] = tidyrbcs(rbc)
% assigns the boundary conditions. Returns them as function handle.

if isdouble(rbc)
    if numel(rbc)==1
        rbcfh = @(u) u(1) - rbc; 
    else
       error('rbc:Input','Left boundary condition is not of robin-type'); 
    end
elseif isa(rbc,'function_handle')
    rbcfh = rbc; % assume that it is in the correct format. 
else
   error('rbc:Input','Left boundary condition is not function handle or double'); 
end

end





