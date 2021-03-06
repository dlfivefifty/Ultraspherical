function u = mldivide(N,rhs,varargin)
% \ mldivide, backslash
%
% MLDIVIDE(N,RHS) solves the linear system N*u=RHS, for u,  where N
% is a ultraop object and rhs is either a scalar or a chebfun
% object. If rhs is a scalar then the rhs is converted into constant valued
% chebfun.

cfs = N.num;
if(N.order==1)
    %either lbc or rbc needs to be set.
    if(isempty(N.lbc) && isempty(N.rbc))
        warning('The operator may not have enough boundary conditions');
    end
elseif(N.order == 2)
    %both lbc and rbc need to be set.
    if( isempty(N.lbc) || isempty(N.rbc))
        warning('The operator may not have enough boundary conditions');
    end
end

if(~isa(rhs,'chebfun'))
    %try and get one.
    crhs = chebfun(rhs);
else
    crhs=rhs;
end

if(isinf(cfs))
    % Adaptive call to solver.
    cfs = max(10,N.order+1);
    rhs = prolong(fliplr(chebpoly(crhs)),cfs);
    rhs = rhs(:);
    
    %map to the right space.
    s=0;
    while(s<N.order)
        %convert rhs to correct coeffs.
        rhs = transMat(cfs,s)*rhs;
        s=s+1;
    end
%     if(N.order==1)
%         if(~isempty(N.lbc))
%             if(iscell(N.lbc))
%                 lvalue = N.lbc{1};
%             else
%                 lvalue = N.lbc;
%             end
%             rhs(end) = lvalue;
%         else
%             if(iscell(N.rbc))
%                 rvalue = N.rbc{1};
%             else
%                 rvalue = N.rbc;
%             end
%             rhs(end) = rvalue;
%         end
%     elseif(N.order==2)
%         % both lbc and rbc must not be empty.
%         if(iscell(N.lbc))
%             lvalue = N.lbc{1};
%         else
%             lvalue = N.lbc;
%         end
%         if(iscell(N.rbc))
%             rvalue = N.rbc{1};
%         else
%             rvalue = N.rbc;
%         end
%         rhs(end-1)=lvalue;
%         rhs(end) = rvalue;
%     elseif(N.order>2)
%         bcvalues = N.lbc;
%         rhs(end-N.order+1:end)=bcvalues;
%     end
    [A bcvalues]=realisation(N,cfs);
    rhs(end-length(bcvalues)+1:end) = bcvalues; 
    
    u = A\rhs;

    %while(max(abs(u(floor(length(u)/1.2):end)))>eps && cfs<20000)
    while(max(abs(u(end-9:end)))>eps && cfs<60000)
        %if(nargin>2&&cfs>1000),break; end;  %threshold to go to gmres.
        cfs=floor(sqrt(2)*cfs);
        rhs = prolong(fliplr(chebpoly(crhs)),cfs);
        rhs = rhs(:);
        
        %map to the right space.
        s=0;
        while(s<N.order)
            %convert rhs to correct coeffs.
            rhs = transMat(cfs,s)*rhs;
            s=s+1;
        end
        rhs(end-length(bcvalues)+1:end) = bcvalues;
%         if(N.order==1)
%             if(~isempty(N.lbc))
%                 if(iscell(N.lbc))
%                     lvalue = N.lbc{1};
%                 else
%                     lvalue = N.lbc;
%                 end
%                 rhs(end) = lvalue;
%             else
%                 if(iscell(N.rbc))
%                     rvalue = N.rbc{1};
%                 else
%                     rvalue = N.rbc;
%                 end
%                 rhs(end) = rvalue;
%             end
%         elseif(N.order==2)
%             % both lbc and rbc must not be empty.
%             if(iscell(N.lbc))
%                 lvalue = N.lbc{1};
%             else
%                 lvalue = N.lbc;
%             end
%             if(iscell(N.rbc))
%                 rvalue = N.rbc{1};
%             else
%                 rvalue = N.rbc;
%             end
%             rhs(end-1)=lvalue;
%             rhs(end) = rvalue;
%         elseif(N.order>2)
%             bcvalues = N.lbc;
%             rhs(end-N.order+1:end)=bcvalues;
%         end
        %rhs = (1./(1:length(rhs)))'.*rhs;
        A = realisation(N,cfs);%*diag(1./(1:cfs));
        %if(cfs<2000), A=full(A);end
        %if(cfs>4000), cfs, end
        %u = (A*diag(1./(1:cfs)))\rhs;
        %u=diag(1./(1:cfs))*u;
        %defaultParms = spparms('autoamd',0);
        %spparms('autommd',0);
        u = A\rhs;
       % spparms('autoamd',0);
        %u = realisation(N,cfs)\rhs;
    end
elseif(isfinite(cfs))
    % Non-adaptive call to solver.
    rhs=crhs;
    
    % form rhs.   
    rhs = prolong(fliplr(chebpoly(rhs)),cfs);
    rhs = rhs(:);
    
    %map to the right space.
    s=0;
    while(s<N.order)
        %convert rhs to correct coeffs.
        rhs = transMat(cfs,s)*rhs;
        s=s+1;
    end
    if(N.order==1)
        if(~isempty(N.lbc))
            rhs(end) = N.lbc;
        else
            rhs(end) = N.rbc;
        end
    elseif(N.order==2)
        % both lbc and rbc must not be empty.
        rhs(end-1)=N.lbc;
        rhs(end) = N.rbc;
    elseif(N.order>2)
        bcvalues = N.lbc;
        rhs(end-N.order+1:end)=bcvalues;
    end
    %Solve system
    A = realisation(N,cfs);%*diag(1./(1:cfs));
    u = A\rhs;
end

%u = truncate(u,eps);

%Scale tiny vector back.
%u = (1./(1:length(u)))'.*u;

u = truncate(u,eps);

%form chebfun as solution.
u = chebfun(flipud(u),'coeffs');
end