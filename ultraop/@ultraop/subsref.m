% SUBSREF,  subassign for an ultraop object. 
% 
% N(n) returned an n by n discretisation of the ultraop N. 
% 
% N.property returns the value of the property for the ultraop N. 


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