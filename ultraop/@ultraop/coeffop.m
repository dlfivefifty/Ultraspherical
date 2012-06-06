% ULTRAOP Constructor for ultraop
% 
% ULTRAOP(F) constructs a ultraop object. 

classdef coeffop
    properties
        order  %Order of DE
        lbc    %Left boundary condition
        rbc    %Right boundary condition
        DEcoeffs %Coefficients of DE (quasi-matrix of chebfuns)
        num   %Solving DE nonadaptively with num points. 
    end
    methods
        function N = coeffop(DE,varargin)
            % This is the coeffop constructor. It sets up a coeffop object
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
                   N.lbc = varargin{k}; 
                   k=k+2;
               elseif(strcmpi('rbc',varargin{k-1}))
                   N.rbc = varargin{k}; 
                   k=k+2;
               else
                   k=k+1;
               end
            end
            
            if(iscell(DE))
                % cell array of function handles
                for j=1:length(DE)
                    DEcheb(:,j) = chebfun(DE{j});
                    % compute the true polynomial interpolant.
                    DEcheb(:,j) = chebfun(DE{j},length(DEcheb(:,j)));
                    %DEcheb(:,j) = chebfun(DE{j},2*length(DEcheb(:,j))); 
                end
                N.DEcoeffs = DEcheb; DE = DEcheb;
            end
            
            if(~isa(DE(:,end),'chebfun'))
                warning('DE should be a cell array or a quasi-matrix of chebfuns'); 
            end
            
            % Construct the coeffop. 
            if(N.order==0)
                if(~isempty(DE))
                    %User wants order to be computed. 
                    N.order = size(DE,2);
                    N.DEcoeffs = DE;
                else
                    %give back empty coeffop. 
                    return;
                end
            end
            if(N.order>0)
                % User wants that order, and default (=0) DEcoeffs. 
                chebzero = chebfun(@(x) 0*x);
                %chebzero = @(x) 0*x;
                while(size(DE,2)<N.order)
                    DE =[chebzero DE];
                end
                if(size(DE,2)>N.order+1)
                    %truncate the end of DE (higher order coeffs)
                    DE = DE(:,1:N.order+1);
                end
            end
            N.DEcoeffs = DE;
        end    
        function u = mldivide(N,rhs,varargin)
            % MLDIVIDE(N,RHS) This solves the linear system N*u=RHS where N
            % is a coeffop object and rhs is either a scalar or a chebfun
            % object. 
            
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
                    if(N.order==1)
                        if(~isempty(N.lbc))
                            if(iscell(N.lbc))
                                lvalue = N.lbc{1};
                            else
                                lvalue = N.lbc;
                            end
                            rhs(end) = lvalue;
                        else
                            if(iscell(N.rbc))
                                rvalue = N.rbc{1};
                            else
                                rvalue = N.rbc;
                            end
                            rhs(end) = rvalue;
                        end          
                    elseif(N.order==2)
                        % both lbc and rbc must not be empty.
                        if(iscell(N.lbc))
                                lvalue = N.lbc{1};
                            else
                                lvalue = N.lbc;
                        end
                        if(iscell(N.rbc))
                            rvalue = N.rbc{1};
                        else
                            rvalue = N.rbc;
                        end
                        rhs(end-1)=lvalue; 
                        rhs(end) = rvalue;
                    elseif(N.order>2)
                        bcvalues = N.lbc;
                        rhs(end-N.order+1:end)=bcvalues;
                    end
                u = realisation(N,cfs)\rhs;
                
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
                    if(N.order==1)
                        if(~isempty(N.lbc))
                            if(iscell(N.lbc))
                                lvalue = N.lbc{1};
                            else
                                lvalue = N.lbc;
                            end
                            rhs(end) = lvalue;
                        else
                            if(iscell(N.rbc))
                                rvalue = N.rbc{1};
                            else
                                rvalue = N.rbc;
                            end
                            rhs(end) = rvalue;
                        end          
                    elseif(N.order==2)
                        % both lbc and rbc must not be empty.
                        if(iscell(N.lbc))
                            lvalue = N.lbc{1};
                        else
                            lvalue = N.lbc;
                        end
                        if(iscell(N.rbc))
                            rvalue = N.rbc{1};
                        else
                            rvalue = N.rbc;
                        end
                        rhs(end-1)=lvalue; 
                        rhs(end) = rvalue;
                    elseif(N.order>2)
                        bcvalues = N.lbc;
                        rhs(end-N.order+1:end)=bcvalues;
                    end
                    %rhs = (1./(1:length(rhs)))'.*rhs;
                    A = realisation(N,cfs);%*diag(1./(1:cfs));
                    if(cfs<2000), A=full(A);end
                    if(cfs>4000), cfs, end
                    %u = (A*diag(1./(1:cfs)))\rhs;
                    %u=diag(1./(1:cfs))*u;
                    u = A\rhs;
                     %u = realisation(N,cfs)\rhs;
                end
            elseif(isfinite(cfs))
                % Non-adaptive call to solver. 
                rhs=crhs;
                %T2U = transMat(2*cfs,0); sT2U = transMat(cfs,0);
                %U2P = transMat(2*cfs,1); sU2P = transMat(cfs,1);
                
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
        function A = realisation(N,cfs)
            % REALISATION(N,CFS) This provides a realisation of the coeffop
            % object with size CFS. It returns a sparse matrix which can be
            % used to solve a collocation problem with CFS collocation
            % pts. 
            
             DE = fliplr(N.DEcoeffs);

                %Form highest order
                s=N.order-1;
                A=speye(cfs);
                while(s>=0)
                    A = A*diffMat(cfs,s);
                    s=s-1;
                end
                if(size(DE,2)>N.order)
                    leading=fliplr(chebpoly(DE(:,end)));
                    if(length(leading)==1)
                        A = leading*A;
                    elseif(length(leading)>1)
                        A = MultMat(leading,cfs,N.order)*A;
                        A = A(1:cfs,:);
                    else
                        error('Only constant leading terms\n');
                    end
                    DE = DE(:,1:end-1);
                end
                
                for j=1:N.order
                    %form D^(j-1) term.
                    cheb = chebpoly(DE(:,j)); cheblen=length(cheb);
                    Cheb = prolong(fliplr(cheb),cfs);
                    Cheb = Cheb(:);
                    %zero coefficient
                    if(cheblen==0),break;end
                    if(cheblen>1)
                        M1 = MultMat(Cheb,cfs,j-1);
                    else
                       M1 = Cheb(1)*speye(cfs,cfs);
                    end
                        s=j-1;
                        while(1+s<=N.order)
                            M1=transMat(cfs,s)*M1;
                            s=s+1;
                        end
                    %truncate.
                    %M1 = M1(1:cfs,:);

                    % Multiply with diffMat
                    s=-1;
                    while(j+s>0)
                        M1 = M1*diffMat(cfs,j+s-1);
                        s=s-1;
                    end
                    A = A + M1;
                    %A =  diag(1./(1:cfs))*(A+M1);
                end
                
                %Set boundary conditions. 
                if(N.order==1)
                    if(~isempty(N.lbc))
                        %LBC matrix row.
                        if(iscell(N.lbc))
                            if(strcmpi(N.lbc{2},'dirichlet'))
                                BCrow = kron(ones(1,cfs),[1 -1]); 
                                BCrow = BCrow(1:cfs);
                                %BCrow = BCrow(1:cfs).*(1./(1:cfs));
                            elseif(strcmpi(N.lbc{2},'neumann'))
                                BCrow = (0:cfs-1).^2.*(-1).^(1:cfs);
                            end
                        else
                            BCrow = kron(ones(1,cfs),[1 -1]); 
                            BCrow = BCrow(1:cfs);
                            %BCrow = BCrow(1:cfs).*(1./(1:cfs));
                        end
                        %wipe one row of A and insert BCrow.
                        A(cfs,:)=BCrow;
                    elseif(~isempty(N.rbc))
                        %RBC matrix row.
                        if(iscell(N.rbc))
                            if(strcmpi(N.rbc{2},'dirichlet'))
                                BCrow = ones(1,cfs);
                            elseif(strcmpi(N.rbc{2},'neumann'))
                                BCrow = (0:cfs-1).^2;
                            end
                        else
                            BCrow = ones(1,cfs);
                        end
                        %wipe one row of A and insert BCrow.
                        A(cfs,:)=BCrow;
                    end

                   
                elseif(N.order==2)
                    % both lbc and rbc must not be empty.
                    if(iscell(N.lbc))      
                        if(strcmpi(N.lbc{2},'dirichlet'))
                            LBC = kron(ones(1,cfs),[1 -1]); 
                            LBC = LBC(1:cfs);
                        elseif(strcmpi(N.lbc{2},'neumann'))
                            LBC = (0:cfs-1).^2.*(-1).^(1:cfs);
                        end
                        A(cfs-1,:)=LBC;
                    elseif(~isempty(N.lbc))
                         LBC = kron(ones(1,cfs),[1 -1]); 
                         LBC = LBC(1:cfs);
                         A(cfs-1,:)=LBC;
                    end
                        %RBC matrix row.
                        if(iscell(N.rbc))
                            if(strcmpi(N.rbc{2},'dirichlet'))
                                RBC = ones(1,cfs);
                            elseif(strcmpi(N.rbc{2},'neumann'))
                                RBC = (0:cfs-1).^2;
                            end
                        elseif(~isempty(N.rbc))
                            RBC = ones(1,cfs);
                            A(cfs,:)=RBC;
                        end 
                    %wipe two rows of A. 
                    
                elseif(N.order>2)
                   %Need to set N.order boundary conditions.
                   %choose them at the chebyshev points?
                   BC=cos(bsxfun(@times,(0:cfs-1),acos(chebpts(N.order))));
                   A(cfs-N.order+1:cfs,:)=BC;
                end
        end
        function y=subsref(N,int)
              % Implement a special subscripted assignment
                switch int(1).type
                case '()'
                    ind = int.subs{:};
                    y = N.realisation(ind);
                case '.'
                    if length(int)>1
                        y = N.(int(1).subs)(int(2).subs{:});
                    else
                        y = N.(int.subs);
                    end
                end
        end
        function y = mtimes(N,y)
           % Computes N*y. 
           if(~isa(y,'chebfun'))
               error('Coeffop:mtimes','Not a chebfun\n');
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
    end %end methods.   
end

