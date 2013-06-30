function [A bcvalues] = realisation(N,cfs,varargin)
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

% %Set boundary conditions.
% if(N.order==1)
%     if(~isempty(N.lbc))
%         %LBC matrix row.
%         if(iscell(N.lbc))
%             if(strcmpi(N.lbc{2},'dirichlet'))
%                 BCrow = kron(ones(1,cfs),[1 -1]);
%                 BCrow = BCrow(1:cfs);
%                 %BCrow = BCrow(1:cfs).*(1./(1:cfs));
%             elseif(strcmpi(N.lbc{2},'neumann'))
%                 BCrow = (0:cfs-1).^2.*(-1).^(1:cfs);
%             end
%         else
%             BCrow = kron(ones(1,cfs),[1 -1]);
%             BCrow = BCrow(1:cfs);
%             %BCrow = BCrow(1:cfs).*(1./(1:cfs));
%         end
%         %wipe one row of A and insert BCrow.
%         A(cfs,:)=BCrow;
%     elseif(~isempty(N.rbc))
%         %RBC matrix row.
%         if(iscell(N.rbc))
%             if(strcmpi(N.rbc{2},'dirichlet'))
%                 BCrow = ones(1,cfs);
%             elseif(strcmpi(N.rbc{2},'neumann'))
%                 BCrow = (0:cfs-1).^2;
%             end
%         else
%             BCrow = ones(1,cfs);
%         end
%         %wipe one row of A and insert BCrow.
%         A(cfs,:)=BCrow;
%     end
%     
%     
% elseif(N.order==2)
%     % both lbc and rbc must not be empty.
%     if(iscell(N.lbc))
%         if(strcmpi(N.lbc{2},'dirichlet'))
%             LBC = kron(ones(1,cfs),[1 -1]);
%             LBC = LBC(1:cfs);
%         elseif(strcmpi(N.lbc{2},'neumann'))
%             LBC = (0:cfs-1).^2.*(-1).^(1:cfs);
%         end
%         A(cfs-1,:)=LBC;
%     elseif(~isempty(N.lbc))
%         LBC = kron(ones(1,cfs),[1 -1]);
%         LBC = LBC(1:cfs);
%         A(cfs-1,:)=LBC;
%     end
%     %RBC matrix row.
%     if(iscell(N.rbc))
%         if(strcmpi(N.rbc{2},'dirichlet'))
%             RBC = ones(1,cfs);
%         elseif(strcmpi(N.rbc{2},'neumann'))
%             RBC = (0:cfs-1).^2;
%         end
%     elseif(~isempty(N.rbc))
%         RBC = ones(1,cfs);
%         A(cfs,:)=RBC;
%     end
%     %wipe two rows of A.
%     
% elseif(N.order>2)
%     %Need to set N.order boundary conditions.
%     %choose them at the chebyshev points?
%     BC=cos(bsxfun(@times,(0:cfs-1),acos(chebpts(N.order))));
%     A(cfs-N.order+1:cfs,:)=BC;
% end

if nargin < 3
    
if isempty(N.lbc) && isempty(N.rbc), bcvalues=[]; return; end; 

bcvalues=[];
rbcrow = cfs;

if ~isempty(N.rbc)
    if isa(N.rbc,'double')
        A(rbcrow,:) = ones(1,cfs); rbcrow=rbcrow-1;
        bcvalues = [bcvalues N.rbc];
    elseif isa(N.lbc,'function_handle')
        [Wa Wb r] = recoverCoeffsBC(N);
        rows = BCrows(Wa,Wb,cfs);
        A(end-length(r)+1:end,:) = rows;
        bcvalues = r;
        return; 
    else
        error('RBC:INPUT','RBC not of correct type');
    end   
end


if ~isempty(N.lbc)
    if isa(N.lbc,'double')
        A(rbcrow,:) = (-1).^(0:cfs-1); 
        bcvalues = [N.lbc bcvalues];
    elseif isa(N.lbc,'function_handle')
        [Wa Wb r] = recoverCoeffsBC(N);
        rows = BCrows(Wa,Wb,cfs);
        A(end-length(r)+1:end,:) = rows;
        bcvalues = r;
        return; 
    else
        error('LBC:INPUT','LBC not of correct type');
    end       
end
        
end

end