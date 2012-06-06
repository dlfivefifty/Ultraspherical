N = coeffop({@(x) x, @(x) cos(x)});
N.lbc = 1; N.rbc =1;
u=N\0;

time = [];
Nint = [];
n=10;
while(n<10000)
    N = coeffop({@(x) x, @(x) cos(x)},n);N.lbc = 1; N.rbc =1;
    s=tic; u=N\0; time = [time toc(s)];
    Nint=[Nint n]; n=2*n;
end
loglog(Nint,time,'r','LineWidth',2); hold on;
loglog(1e3:1e5,1e-6*(1e3:1e5).^1.5,'k--','LineWidth',2);
%%
ctime = []; dtime=[];
cNint = [];

n=10;
while(n<1000)
    N = chebop(@(x,u) diff(u,2) + x.*diff(u,1) + cos(x).*u);
    N.lbc = 1; N.rbc =1;
    %s=tic; u=N(n)\zeros(n,1); ctime = [ctime toc(s)];
    S=svd(N(n));
    ctime = [ctime 1./S(end)];
    dtime = [dtime S(1)];
    cNint=[cNint n]; n=2*n;
end


loglog(cNint,ctime,'b','LineWidth',2); hold on;
loglog(cNint,dtime,'r','LineWidth',2); hold on;
%loglog(5e2:5e3,1e-9*(5e2:5e3).^3,'k--','LineWidth',2);

title('Condition of a dense spectral matrix','FontSize',16);
%legend('Coeffop','O(n^{1.5})','Chebfun','O(n^3)','Location','best');
xlabel('size'); %ylabel('time');
legend('||A||','||A^{-1}||','Location','best')