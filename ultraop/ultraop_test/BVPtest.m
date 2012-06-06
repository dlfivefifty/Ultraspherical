function Passed = BVPtest()
  % Solve 2nd order BVP u''(x) + a(x)u'(x) + u(x) = f(x), with bc
  
  % Solve u''(x) + sin(x)u'(x) + cos(x)u(x) = 0 with u'(-1)=1, u(1) =1.
  N = chebop(@(x,u) diff(u,2)+sin(x).*diff(u) + cos(x).*u); 
  N.lbc=@(u) diff(u)-1; N.rbc=1;
  u = N\0; 
  N = ultraop({@(x) 1+0*x,@(x) sin(x),@(x) cos(x)}); 
  N.lbc={1,'neumann'}; N.rbc=1;
  y = N\0; 
  Pass(1) = (norm(u-y)<2e-11);

  % Solve u''(x) + cos(x)u(x) = 0 with u(-1)=1, u(1) =1, compare different
  % syntax.
  N = ultraop({@(x) 1+0*x,@(x) 0+0*x,@(x) cos(x)}); 
  N.lbc=1; N.rbc=1;
  u = N\0; 
  N = ultraop({@(x) 1+0*x,@(x) 0+0*x,@(x) cos(x)},'order',2); 
  N.lbc=1; N.rbc=1;
  y = N\0; 
  Pass(2) = (norm(u-y)==0);
  
  % Coefficients given as chebfuns 
  x = chebfun('x'); 
%   a = 0*x; b=chebfun(@(x) cos(x));
  N = ultraop({1+0*x,0*x,cos(x)}); 
  N.lbc=1; N.rbc=1;
  u = N\0;  
  Pass(3) = (norm(u-y)==0);

  if(all(Pass))
     Passed=1; 
  else
     Passed=0;
  end
  
  %% to do
  % finish display. 
  % error check when N.num is finite! ! ! !
  % Do order=4 stuff. 
  
% N = coeffop({@(x)  5e6*x.*(x.^2-0.5),@(x)  15e6*(x.^2-0.5)});
% N.lbc=-2;
% N.rbc=4;
% y=N\0;
end