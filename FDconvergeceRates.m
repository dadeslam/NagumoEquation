%% Analysis of the right convergence rates in space and time of the schemes

clear all
close all
clc

%The equation to be solved is u_t=u_xx+u(1-u)(a-u) with a in [0,1/2).
%and some suitable boundary condition.

% First of all let's compare 2 different schemes: Pseudo Crank-Nicolson and Explicit scheme
% Making the analyisis of the convergence rate and time required to attain a
% certain tolerance using the exact solution.

a = 0.25;

N = 1001; %Space discretization points

L = 100;
T = 0.1;

x = linspace(-L,L,N)';
dx = x(2)-x(1);

c = 0.2; %So for stability we set dt = 2*c*dx^2/(4+a*dx^2)
dt = min(T/100,c*dx^2/(4+a*dx^2));
time = (0:dt:T);
M = length(time);

r = dt/(2*dx^2);

EM = diag((1-4*r)*ones(N,1),0) + diag(2*r*ones(N-1,1),1) + diag(2*r*ones(N-1,1),-1); %Matrix for the explicit scheme:
% u_new = EM*u_now +dt*f(u_now)
EM(1,1:2) = [1,0];
EM(N,N-1:N) = [0,1];

a1 = 2*r*ones(N,1);
a2 = -r*ones(N-1,1);
A1 = diag(a1);
A2 = diag(a2,1);
CN = A1+A2+A2'; %Matrix required for the Pseudo Crank-Nicolson method
%(I+CN)u_n+1 = (I-CN)*u_n+dt*f(u_n)


%Reference exact solution:
A = 1;
B= 1;
g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
uexact = (A*exp(g1)+a*B*exp(g2))./(A*exp(g1)+B*exp(g2)+1); %Exact solution at time T

%Boundary conditions for this test are u(xL,t) = 0 and u(xR,t) = 1. 

    
% Initial condition for this test
u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
uE = u0; %initial solution for the explicit scheme
uC = u0; %initial solution for the Pseudo Crank-Nicolson scheme

B = (speye(N)+CN); %Symmetric tridiagonal matrix used in the scheme
m = 10^15;
B(1,1) = m;
B(end,end) = m;
R = chol(B);

for t = time(2:end)
    b = (eye(N)-CN)*uC+dt*uC.*(1-uC).*(uC-a);
    %Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1
    b(1) = 0;
    b(end) = m;
    uC = R'\b;
    uC = R\uC;

end

for t = time(2:end)
    uE(1) = 0;
    uE(end) = 1;
    uE = EM*uE++dt*uE.*(1-uE).*(uE-a);
end


figure;
plot(x,uE,'r*',x,uC,'ko',x,uexact,'g-')
title('Plot of the 3 solutions to get a first idea')
legend('Explicit solution','ps. Crank-Nicolson solution','Exact solution')


%% Space convergence rate analysis
L = 100;
T = 0.1;
a = 0.25;

dt = T/1001;
time = (0:dt:T);
M = length(time);

xrange = (500:100:900);
count = 1;
errorE = zeros(length(xrange),1);

for N = xrange
    x = linspace(-L,L,N)';
    dx = x(2)-x(1);


    r = dt/(2*dx^2);

    EM = diag((1-4*r)*ones(N,1),0) + diag(2*r*ones(N-1,1),1) + diag(2*r*ones(N-1,1),-1); %Matrix for the explicit scheme:
    % u_new = EM*u_now +dt*f(u_now)
    EM(1,1:2) = [1,0];
    EM(N,N-1:N) = [0,1];


    %Reference exact solution:
    A = 1;
    B= 1;
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (A*exp(g1)+a*B*exp(g2))./(A*exp(g1)+B*exp(g2)+1); %Exact solution at time T

    %Boundary conditions for this test are u(xL,t) = 0 and u(xR,t) = 1. 


    % Initial condition for this test
    u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
    uE = u0; %initial solution for the explicit scheme

    for t = time(2:end)
        uE(1) = 0;
        uE(end) = 1;
        uE = EM*uE++dt*uE.*(1-uE).*(uE-a);
    end
    errorE(count) = norm(uE-uexact,inf);
    count = count + 1;
end

label = 'Tested number of space points with dt='+string(dt);

ord1 = (xrange/xrange(1)).^(-1)*errorE(1);
ord2 = (xrange/xrange(1)).^(-2)*errorE(1);
ord3 = (xrange/xrange(1)).^(-3)*errorE(1);

figure;
loglog(xrange,errorE,'r*',xrange,ord1,'b-',xrange,ord2,'k-',xrange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in space measured at t=0.1 for Explicit scheme');
xlabel(label);
ylabel('Maximum error');

%Space analysis for pseudo CN
dt = T/401;
count = 1;
time = (0:dt:T);
errorC = zeros(length(xrange),1);
for N = xrange
    
    dx = 2*L/(N-1);    
    x = (-L:dx:L)';
    A = 1;
    B = 1;
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (A*exp(g1)+a*B*exp(g2))./(A*exp(g1)+B*exp(g2)+1); %Exact solution at time T
    
    M = length(time);

    r = dt/(2*dx^2);
    a1 = 2*r*ones(N,1);
    a2 = -r*ones(N-1,1);
    A = diag(a1) + diag(a2,1) + diag(a2,-1);

    %The scheme writes, in matrix form, as follows:
    %(I+A)u_n+1 = (I-A)*u+dt*u.*(1-u).*(u-a)

    u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
    uC = u0;
    m = 10^15;
    B = (speye(N)+A);
    B(1,1) = m;
    B(end,end) = m;
    R = chol(B);

    for t = time(2:end)
        b = (eye(N)-A)*uC+dt*uC.*(1-uC).*(uC-a);
        b(1) = 0;
        b(end) = m;

        uC = R'\b;
        uC = R\uC;
    end
    
    errorC(count) = norm(uC-uexact,inf);
    count = count + 1;
end


label = 'Tested number of space points with dt='+string(dt);

ord1 = (xrange/xrange(1)).^(-1)*errorC(1);
ord2 = (xrange/xrange(1)).^(-2)*errorC(1);
ord3 = (xrange/xrange(1)).^(-3)*errorC(1);

figure;
loglog(xrange,errorC,'r*',xrange,ord1,'b-',xrange,ord2,'k-',xrange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in space measured at t=0.1 for ps. Crank-Nicolson');
xlabel(label);
ylabel('Maximum error');

%% Time convergence rate analysis
L = 100;

T = 0.1;
a = 0.25;

trange = T./(50:10:90);

errorE = zeros(length(trange),1);
dx = sqrt(4*trange(1)/(2-a*trange(1)));


x = (-L:dx:L)';

count = 1;

%Explicit scheme
for dt = trange
    
    time = (0:dt:T);
    A = 1;
    B= 1;
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (A*exp(g1)+a*B*exp(g2))./(A*exp(g1)+B*exp(g2)+1); %Exact solution at time T
    
    N = length(x);
    M = length(time);
    r = dt/(2*dx^2);
        
    %Matrix for the explicit scheme
    d1 = (1-4*r)*ones(N,1);
    d2 = 2*r*ones(N-1,1);
    EM = diag(d1) + diag(d2,1) + diag(d2,-1);

    u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
    uE = u0;
    
    for t = time(2:end) 
        uE(1) = 0;
        uE(end) = 1;
        uE = EM*uE+dt*uE.*(1-uE).*(uE-a);
    end
    
    errorE(count) = norm(uE-uexact,inf);
    count = count + 1;
end

trange = T./trange;

ord1 = (trange/trange(1)).^(-1)*errorE(1);
ord2 = (trange/trange(1)).^(-2)*errorE(1);
ord3 = (trange/trange(1)).^(-3)*errorE(1);
 
label = 'Tested number of points in time for dx='+string(dx);

figure;
loglog(trange,errorE,'r*',trange,ord1,'b-',trange,ord2,'k-',trange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in time measured at t=0.1 for Explicit scheme');
xlabel(label);
ylabel('Maximum error');


%Pseudo Crank-Nicolson
N = 3000;

x = linspace(-L,L,N)';
dx = x(2)-x(1);

c = 0.2; %So for stability we set dt = 2*c*dx^2/(4+a*dx^2)

trange = 15:5:30;
errorC = zeros(length(trange),1);
count = 1;

for M = trange
    dt = T/(M-1);
    time = (0:dt:T);

    r = dt/(2*dx^2);
    a1 = 2*r*ones(N,1);
    a2 = -r*ones(N-1,1);
    A1 = diag(a1);
    A2 = diag(a2,1);
    CN = A1+A2+A2'; %Matrix required for the Pseudo Crank-Nicolson method
    %(I+CN)u_n+1 = (I-CN)*u_n+dt*f(u_n)


    %Reference exact solution:
    A = 1;
    B= 1;
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (A*exp(g1)+a*B*exp(g2))./(A*exp(g1)+B*exp(g2)+1); %Exact solution at time T

    %Boundary conditions for this test are u(xL,t) = 0 and u(xR,t) = 1. 


    % Initial condition for this test
    u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
    uC = u0; %initial solution for the Pseudo Crank-Nicolson scheme

    B = (speye(N)+CN); %Symmetric tridiagonal matrix used in the scheme
    m = 10^15;
    B(1,1) = m;
    B(end,end) = m;
    R = chol(B);

    for t = time(2:end)
        b = (eye(N)-CN)*uC+dt*uC.*(1-uC).*(uC-a);
        %Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1
        b(1) = 0;
        b(end) = m;
        uC = R'\b;
        uC = R\uC;

    end
    errorC(count) = norm(uC-uexact,inf);
    count = count + 1;
end

ord1 = (trange/trange(1)).^(-1)*errorC(1);
ord2 = (trange/trange(1)).^(-2)*errorC(1);
ord3 = (trange/trange(1)).^(-3)*errorC(1);

label = 'Tested number of points in time for dx='+string(dx);

figure;
loglog(trange,errorC,'r*',trange,ord1,'b-',trange,ord2,'k-',trange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in time measured at t=0.1 for ps. Crank-Nicolson');
xlabel(label);
ylabel('Maximum error');

% % %So up to now we have just seen that from the point of view of
% % %computational speed, in terms of accuracy they are equivalent and in terms
% % %of stability we should prefer pseudo CN, because unconditionally stable. So we
% % %will make the following analysis with pseudo CN method.

%% Just for a sake of completeness here it is the implementation of a numerical scheme
%% which is second order both in space and time, it is not a finite differences scheme
%% but after a finite difference discretization in space (of second order) we solve the
%% pseudo discretized problem with a second order time integration, for instance a Runge Kutta 2.
%% So what is actually implemented is METHOD OF LINES

a = 0.25;
T = 0.1;   
Mrange = 2.^(6:9);
N = 6002;
H = 0.1; 
b = 0;
L = 100;
x = linspace(-L,L,N)';
dx = x(2)-x(1);

%Reference exact solution:
A = 1;
B= 1;
g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
uexact = (A*exp(g1)+a*B*exp(g2))./(A*exp(g1)+B*exp(g2)+1); %Exact solution at time T

%Boundary conditions for this test are u(xL,t) = 0 and u(xR,t) = 1. 
    
% Initial condition for this test
u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);

D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/dx^2,1,N));

% So that u_xx = D2 * u
% Now instead of u(x,t) we have u=u(t) = [u1(t),...,uN(t)]
% and the equation, just in time, becomes u'(t) =
% D2*u(t)+u(t)*(1-u(t))*(u(t)-a)

f = @(u) D2*u+u.*(1-u).*(u-a);

error = zeros(length(Mrange),1);
count = 1;

for M = Mrange
    time = linspace(0,T,M)';
    u = u0;
    dt = time(2)-time(1);
    for t = time(2:end)
        u = u + dt/2*(f(u)+f(u + dt*f(u)));
        u(1) = 0;
        u(end) = 1;
    end
    error(count) = norm(u-uexact,inf);
    count = count + 1;
end

ord1 = (Mrange/Mrange(1)).^(-1)*error(1);
ord2 = (Mrange/Mrange(1)).^(-2)*error(1);
ord3 = (Mrange/Mrange(1)).^(-3)*error(1);

figure;
loglog(Mrange,error,'r*',Mrange,ord1,'b-',Mrange,ord2,'k-',Mrange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in time measured at t=0.1 for improved Euler method');
xlabel('Tested dt for dx=1e-4');
ylabel('Maximum error');
