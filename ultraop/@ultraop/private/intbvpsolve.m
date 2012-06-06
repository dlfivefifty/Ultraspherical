% function v = intbvpsolve(L,f)
% % Solve a linop equation using integral reformulation
% P = recoverCoeffs(L);             % Obtain the coeff functions of the op
% [Wa Wb r W] = recoverCoeffsBC(L); % Recover the BC coeffs
% v = intbvp(P,f,Wa,Wb,r,W);        % Solve
% 
% end

% function [u,uder] = intbvp(P,f,Wa,Wb,r,W)
% % From Driscoll 2010
% m = size(P,2)-1; % order of ODE
% [C,T,AT,WT,ACm,WCm,H] = intdata(P,Wa,Wb,W);
% % Construct chebfun solution for v = D^m * u, and integration constants.
% v = (ACm - H) \ (f - AT*(WT\r));
% k = WT \ (r-WCm * v);
% % Reconstruct all derivatives of the solution accurately.
% Cjv = v; uder = v;
% for j = 1:m
%     Cjv = C * Cjv;
%     uder(:,j + 1) = Cjv + T(:,1:j)*k(m-j + 1:m);
% end
% u = uder(:,end); % same as u = C^m * v + T * k
% end


% function [C,T,AT,WT,ACm,WCm,H] = intdata(P,Wa,Wb,W)
% % From Driscoll 2010
% m = size(P,2)-1;
% % order of ODE
% d = domain(P(:,1));
% if ~isempty(W), d = union(domain(W),d); end
% a = d(1); b = d(end);
% C = cumsum(d);
% % Formulate quasimatrix for null space of D^m.
% T = chebfun(1,d);
% for k = 2:m, T(:,k) = C * T(:,k-1); end
% % Interior operator A applied to T, by Horner’s method.
% AT = 0; DT = T;
% for j = 1:m
%     AT = AT + diag(P(:,j))*DT;
%     DT = [chebfun(0,d), DT(:,1:m-1)];
% end
% % Boundary operator W applied to T.
% WT = Wa + Wb * toeplitz( [1 zeros(1,m-1)], T(b,:) );
% % Operators applied to C^m integration operator, again using Horner.
% ACm = diag(P (:,1));
% Eb = feval(d,b);
% WCm = Wb(:,1)*Eb;
% for j = 1:m-1
%     ACm = ACm * C + diag(P (:,j + 1));
%     WCm = WCm * C + Wb(:,j + 1)*Eb;
% end
% ACm = ACm * C + diag(P(:,m + 1));
% WCm = WCm * C;
% 
% % Add contributions from non-standard terms
% if ~isempty(W)
%     WCm = [WCm ; W*cumsum(d,m)];
%     WT = [WT ; W*T];
% end
% 
% % Chebop form of outer products.
% op = @(v) AT * (WT\(WCm * v)); % functional expression
% H = linop(@mat,op,d);
% 
%     function m = mat(n)
%         nin = n;
%         [n map breaks] = tidyInputs(n,d,'initdata');      
% 
%         if isempty(map) && isempty(breaks)
%             % Standard case
%             m = AT(chebpts(n,d),:) * (WT\WCm(n)); % n-by-n realization
%         elseif isempty(map)
%             % Breakpoints / no maps 
%             WCmN = feval(get(WCm,'varmat'),nin);
%             m = AT(chebpts(n,d),:) * (WT\WCmN); % n-by-n realization
%         else
%             error('maps are not supported here.')
%         end    
%     end
%     
% end



function [p varargout] = recoverCoeffs(L)
%RECOVERCOEFFS  Recover coefficient functions of a linear operator
% P = RECOVERCOEFFS(L) returns, for a linear operator L, a chebfun
% quasimatrix P such that
%         Lu = P(:,1)*u + P(:,2)*u' + P(:,3)*u" + ... P(:,M+1)*u^(M),
% where M is the difforder of the operator. If L is not linear, an error is
% thrown.
%
% For a block operator L, i.e., one defining a system of equations
%         Lu = [L_{1,1} L_{1,2} ... L_{1,S}] [ u_1 ]
%              [L_{2,1} L_{1,2} ... L_{1,S}] [ u_2 ]
%              [  ...     ...   ...   ...  ] [ ... ]
%              [L_{R,1} L_{R,2} ... L_{R,S}] [ u_S ],
% P will be the RxS cell array such that P{J,K} = RECOVERCOEFFS(L_{J,K}).
%
% [P L] = RECOVERCOEFFS(L) returns also the linop L, which can be useful if
% the input was a linear chebop.
%
% Example 1:
%  [L x] = chebop(@(x,u) 0.5*diff(u,2) - sin(x).*diff(u) + x.*u);
%  p = recoverCoeffs(L)
%  norm(p - [x -sin(x) 0.5])  
%
% Example 2:
%  [L x] = chebop(@(x,u) diff(sin(x).*(diff(cos(x).*u))),[-pi pi]);
%  p = recoverCoeffs(L)
%  norm(p - [-sin(2*x) 1-3*sin(x).^2 sin(2*x)/2])
%
% Example 3:
%  L = chebop(@(x,u,v) [diff(u,2), 0.5*diff(v)+exp(x)]);
%  p = recoverCoeffs(L) 
%  norm([p{:}] - [0 0 1 0 0 0 .5])

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% Convert to linop if input is a chebop. (But don't overwrite input as it's
% more efficient to evaluate the chebop .op than the linearised .oparray!)
if isa(L,'chebop'), L2 = linop(L); else L2 = L; end

% Initialise
s = L2.blocksize;                % Determine the size of the system
m = L2.difforder;                %  and the difforder
x = chebfun('x',L2.fundomain,2); % Construct linear function on the domain
x0 = chebfun(0,L2.fundomain);    %  and the zero function
p = cell(s);                     % Initialise output
p0 = L*repmat(x0,1,s(2));        % Compute non-autonomous component

% The main routine
for hh = 1:s(2)                 % Loop over each of the dependant variables
    x0l = repmat(x0,1,hh-1);    % Set dep vars to the left to zero
    x0r = repmat(x0,1,s(2)-hh); % Set dep vars to the right to zero
    p1 = L*[x0l 1+0*x x0r];     % Evaluate all equations for [0 ... 1 ...0]
    p1 = p1 - p0;               % Subtract non-autonomous compnent
    for ll = 1:s(1)             % Loop over equations and assign
        p{ll,hh} = p1(:,ll);
    end
    xk = x;                            % Update indep var to x
    for kk = 1:max(m(:,hh))            % Loop over each x^k
        tmp = L*[x0l xk(:,kk) x0r]-p0; % Evaluate for u = [0 ... x^k ... 0]
        for ll = 1:s(1)                % Loop over each equation
            if kk > m(ll,hh), continue, end % No coeffs of this order here
            p{ll,hh}(:,kk+1) = tmp(:,ll);   % Assign the ll-th equation
            for jj = 1:kk              % Extract the lower-order terms
                p{ll,hh}(:,kk+1) = p{ll,hh}(:,kk+1) - p{ll,hh}(:,kk+1-jj).*xk(:,jj);
                p{ll,hh}(:,kk+1) = simplify(p{ll,hh}(:,kk+1)); % Simplify
            end
        end
        xk = [xk x.*xk(:,end)/(kk+1)]; % Update indep var to x^k/k!
    end
end

% Tidy the output
if max(s) == 1, p = p{:}; end          % Output quasimatrix if not a system
if nargout == 2, varargout{1} = L; end % Output the linop if required
end

function [Wa Wb r W] = recoverCoeffsBC(L)
%RECOVERCOEFFSBC  Recover coefficient functions of a linear BCs
% [WA WB R] = RECOVERCOEFFS(L) returns, for a linear operator L on a domain 
% [a,b], two matrices WA and WB such that
%   sum_{k=0}^{n} ( WA(:,k)*u^(k)(a) + WB(:,k)*u^(k)(b) ) = R.
%
% [WA WB R W] = RECOVERCOEFFS(L) allows for cases where L has nonstandard 
% or interiooir BCs in its .bc field, in which case W is a row-linop and
%   sum_{k=0}^{n} ( WA(:,k)*u^(k)(a) + WB(:,k)*u^(k)(b) ) = R(1:size(WA,1))
% and
%   W*u = R(size(WA,1)+1:end)
%
% Example 1:
%   N = chebop(@(x,u) diff(u,2),[-1 1]);
%   N.lbc = @(u) u-2;
%   N.rbc = @(u) 3*diff(u)-4;
%   [Wa Wb r] = recoverCoeffsBC(N)
%   Wa =                Wb =                r = 
%          0     1             0     0            2
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

W = [];
if any(~valid)
    idx = find(~valid);
    BC = L.bc;
    
    for type = {'left','right','other'};
        list = BC.(type{:});
        numbc = numel(list);
        while ~isempty(idx) && any(idx(1) <= numbc)
            W = [W ; list(idx(1)).op];
            idx(1) = [];
        end
        idx = idx - numbc;
    end
    r = [r ; c(~valid)];
end

end
    
% function [n map breaks numints] = tidyInputs(n,d,filename)
% % TIDYINPUTS - tidy the inputs to varmat constructors
% %      [n map breaks numints] = tidyInputs(n,d,filename)
% % Ensures that breakpoints are inherited/disregarded correctly and that the
% % lenght of n is right for the given domain/breaks. 
% 
% breaks = []; map = [];
% if iscell(n)
%     if numel(n) > 1, map = n{2}; end
%     if numel(n) > 2, breaks = n{3}; end
%     if isa(breaks,'domain'), breaks = breaks.ends; end
%     n = n{1};
% end
% 
% % Inherit the breakpoints from the domain.
% breaks = union(breaks, d.ends);
% % Throw away breaks (and corresponding n) outside the domain.
% maskr = breaks > d.endsandbreaks(end);     maskl = breaks < d.endsandbreaks(1);       
% if numel(n) > 1,  n(maskr) = [];  n(maskl) = []; end
% breaks(maskl|maskr) = [];
% numints = numel(breaks)-1;
% 
% % Force a default map for unbounded domains.
% if any(isinf(breaks)) && isempty(map)
%     map = maps(domain(breaks));
% end  
% 
% % Tidy up breaks and n.
% if numints == 1
%     % Breaks are the same as the domain ends. Set to [] to simplify.
%     breaks = [];
% elseif numel(breaks) > 2
%     % Repmat n if necessary - or check for mismatch
%     if numel(n) == 1, n = repmat(n,1,numints); end
%     if numel(n) ~= numints
%         if nargin < 3, filename = 'unknown'; end
%         error(['DOMAIN:',filename,':numints'],...
%             'Vector N does not match domain D.');
%     end
% end
% end
    
    






    



