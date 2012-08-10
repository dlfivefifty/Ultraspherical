%% Lee and Greengard's ODE examples

%%
% We try out a selection of ODE examples, using US methods. We do not do
% subdivision of the interval so [1] can do stiffer problems than US
% methods without splitting. However, the US method provides very accurate
% solutions. 

%% EXAMPLE 1
% The first example is
% $$ \varepsilon u''(x) + 2 x u'(x) = 0, \quad u(-1) = 1, ~ u(1) = 1. $$
FS = 'fontsize'; LW = 'linewidth';
ep=1e-4;
N = ultraop(@(x,u) ep*diff(u,2) + 2*x.*diff(u)); N.lbc=-1; N.rbc=1; 
tic, u = N\0; t = toc;
plot(u,'m',LW,2);

%% EXAMPLE 2
% This example is an Airy equation,
% $$ \varepsilon u''(x) - x u(x) = 0, \quad u(-1) = 1, ~ u(1) = 1. $$
ep= 1e-7;
N = ultraop(@(x,u) ep*diff(u,2)-x.*u); N.lbc=1; N.rbc=1;
u = N\0; 
plot(u,'r',LW,2)

%% EXAMPLE 4
% The fourth example is
% $$ \varepsilon u''(x) + (x^2-0.25)u(x) = 0, \quad u(-1) = 1, ~ u(1) = 2. $$
ep = 1e-7; 
N = ultraop(@(x,u) ep*diff(u,2)+(x.^2-0.25).*u); N.lbc=1; N.rbc=2;
u = N\0; 
plot(u,LW,2)

%% EXAMPLE 5. 
% The fifth example is
% $$ \varepsilon u''(x) + x u'(x) - 0.5u(x) = 0, \quad u(-1) = 1, ~ u(1) = 2. $$
ep=1e-5;
N = ultraop(@(x,u) ep*diff(u,2)+x.*diff(u)-0.5*u); N.lbc=1; N.rbc=2;
u = N\0; 
plot(u,LW,2)

%% EXAMPLE 6
% Last example is
% $$ \varepsilon u''(x) - x u'(x) + u(x) = 0, \quad u(-1) = 1, ~ u(1) = 2. $$
ep=1e-4;
N = ultraop(@(x,u) ep*diff(u,2)-x.*diff(u)+u); N.lbc = 1; N.rbc=2; 
u = N\0; 
plot(u,LW,2)

%%
% REFERENCE
%
% [1] J.-Y. Lee and L. Greengard, "A fast adaptive numerical method for
% stiff two-point boundary value problems", SIAM Journal on Scientific
% Computing 18 (1997), 403-429.
%
% [2] Lloyd N. Trefethen, "Lee and Greengard ODE Examples", Chebfun
% example, http://www2.maths.ox.ac.uk/chebfun/examples/ode/, (2012). 
%
