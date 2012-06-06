function Passed = IVPtest()
    % Solve 1st order systems. u'(x) + a(x)u(x) = f(x); with bc.
    
    % u'(x) + u = 0, u(0)=1;
    exact = chebfun(@(x) exp(-1).*exp(-x));   
    N = ultraop({@(x) 1 + 0*x,@(x) 1 + 0*x}); 
    N.lbc=1;
    sol=N\0;
    Pass(1) = (norm(sol-exact)<1e-14);
    
    % u'(x) + 1i*u = 0, u(0)=1; 
    N = chebop(@(x,u) diff(u) + 1i*u); N.lbc=1;
    exact = N\0;
    N = ultraop({@(x) 1 + 0*x,@(x) 1i+0*x}); 
    N.lbc = 1; sol=N\0;
    Pass(2) = norm(exact-sol)<1e-14;
     
    % u'(x) + x*u = 0, u(0)=1;
    N = chebop(@(x,u) diff(u) + x.*u); N.lbc=1;
    exact = N\0;
    N = ultraop({@(x) 1 + 0*x,@(x) x},200); N.lbc=1;
    sol=N\0;
    Pass(3) = (norm(sol-exact)<1e-14);

    % a(x) = 1i*x.
    N = chebop(@(x,u) diff(u) + 1i*x.*u,[-1,1]); N.rbc=1;
    exact = N\0;
    N = ultraop({@(x) 1 + 0*x,@(x) 1i*x}); N.rbc=1;
    sol = N\0;
    Pass(4) = (norm(sol-exact)<1e-14);
    
    % u'(x) + a(x)u(x) = f(x) with bc. 
    N = chebop(@(x,u) diff(u) + cos(x).*u); N.rbc=1;
    f=chebfun(@(x) sin(x)); 
    exact = N\f;
    N = ultraop({@(x) 1 + 0*x,@(x) cos(x)}); N.rbc=1;
    sol = N\f;
    Pass(5) = (norm(sol-exact)<1e-14);
    
    % u'(x) + a(x)u(x) = f(x) with neumann bc. 
    N = chebop(@(x,u) diff(u) + cos(x).*u); N.rbc=@(u) diff(u)-1;
    f=chebfun(@(x) sin(x)); 
    exact = N\f;
    N = ultraop({@(x) 1 + 0*x,@(x) cos(x)}); N.rbc={1,'neumann'};
    sol = N\f;
    Pass(6) = (norm(sol-exact)<1e-12);
    
    if(all(Pass)) 
        Passed=1;
    else
        Passed=0;
    end
end