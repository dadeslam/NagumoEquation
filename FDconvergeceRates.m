%% Analysis of the right convergence rates in space and time of the schemes

clear all
close all
clc

%% Running the code will show all the plots in the associated section, so all the error
%% rates computed analytically.

%% The equation to be solved is u_t=u_xx+u(1-u)(a-u) with a in [0,1/2).
%% and some suitable boundary condition.

% First of all let's compare 2 different schemes: Pseudo Crank-Nicolson and Explicit scheme
% Making the analyisis of the convergence rate

% Then we conclude with the Heun method with centered finite differences

%% Initial plot of the 2 overlapping numerical solutions with the exact one at time T

a = 0.25; %Threshold parameter

N = 1001; %Space discretization points

L = 100; %Length of the positive part of the space domain
T = 0.1; %Final time

x = linspace(-L,L,N)'; %Discretized space domain
dx = x(2)-x(1); %Space discretization interval

c = 0.2; %So for stability we set dt = 2*c*dx^2/(4+a*dx^2), it is important c<=1
dt = min(T/100,c*dx^2/(4+a*dx^2));
time = (0:dt:T); %Discretized time interval
M = length(time);

r = dt/(2*dx^2);

value = dt*a; %Is less than 2 so even Pseudo Crank-Nicolson is stable 

EM = diag((1-4*r)*ones(N,1),0) + diag(2*r*ones(N-1,1),1) + diag(2*r*ones(N-1,1),-1); %Matrix for the explicit scheme:
% u_new = EM*u_now +dt*f(u_now)
EM(1,1:2) = [1,0]; %I change these entries to fix the right BCs
EM(N,N-1:N) = [0,1];

a1 = 2*r*ones(N,1);
a2 = -r*ones(N-1,1);
A1 = diag(a1);
A2 = diag(a2,1);
CN = A1+A2+A2'; %Matrix required for the Pseudo Crank-Nicolson method
%(I+CN)u_n+1 = (I-CN)*u_n+dt*f(u_n)


%Reference exact solution at time t=T:
g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T

%Boundary conditions for this test are u(xL,t) = 0 and u(xR,t) = 1. 

    
% Initial condition for this test
u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
uE = u0; %initial solution for the explicit scheme
uC = u0; %initial solution for the Pseudo Crank-Nicolson scheme

B = (speye(N)+CN); %Symmetric tridiagonal matrix used in the pseudo CN scheme
m = 10^15; %Penalization parameter
B(1,1) = m; %Change the entries to preserve the SPD character but having the right BCs
B(end,end) = m;
R = chol(B); %Cholesky factorization

%Pseudo Crank - Nicolson
for t = time(2:end)
    b = (eye(N)-CN)*uC+dt*uC.*(1-uC).*(uC-a); %Known right hand side for pseudo CN
    %Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1
    b(1) = 0;
    b(end) = m;
    uC = R'\b;
    uC = R\uC;
end

%Explicit Euler with centered finite differences
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

xrange = (500:100:900); %Range of tested number of space points
count = 1;
errorE = zeros(length(xrange),1); %Vector where we will store the errors in space
%Space convergence rate for Explicit Euler with centered finite difference
for N = xrange
    x = linspace(-L,L,N)';
    dx = x(2)-x(1);
    r = dt/(2*dx^2);

    EM = diag((1-4*r)*ones(N,1),0) + diag(2*r*ones(N-1,1),1) + diag(2*r*ones(N-1,1),-1); %Matrix for the explicit scheme:
    % u_new = EM*u_now +dt*f(u_now)
    EM(1,1:2) = [1,0];
    EM(N,N-1:N) = [0,1];


    %Reference exact solution:
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T

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

const = errorE'./xrange.^2; %It gives the constant in time such that ||uNumerical - uExact||<= Const * Delta x^2.
disp('So the average constant for the O(Delta x^2) for Explicit Scheme is:');
c = mean(const)

label = 'Tested number of space points with dt='+string(dt);

ord1 = (xrange/xrange(1)).^(-1)*errorE(1); %Convergence orders to compare with our
ord2 = (xrange/xrange(1)).^(-2)*errorE(1);
ord3 = (xrange/xrange(1)).^(-3)*errorE(1);

figure;
loglog(xrange,errorE,'r*',xrange,ord1,'b-',xrange,ord2,'k-',xrange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in space measured at t=0.1 for Explicit scheme');
xlabel(label);
ylabel('Maximum error');

%Space analysis for pseudo CN
dt = T/401; %We fix a different time interval since here we are more free with the step
count = 1;
time = (0:dt:T);
errorC = zeros(length(xrange),1);
for N = xrange
    
    dx = 2*L/(N-1);    
    x = (-L:dx:L)';
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T
    
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
    R = chol(B); %Usual cholesky factorization

    for t = time(2:end)
        b = (eye(N)-A)*uC+dt*uC.*(1-uC).*(uC-a);
        b(1) = 0;
        b(end) = m;

        uC = R'\b;
        uC = R\uC;
    end
    
    errorC(count) = norm(uC-uexact,inf); %Error for pseudo CN
    count = count + 1;
    
end
value = a*dt; %Is less than 2 so no problem with the stability of Pseudo Crank Nicolson here

const = errorC'./xrange.^2; %It gives the constant in time such that ||uNumerical - uExact||<= Const * Delta x^2.
disp('So the average constant for the O(Delta x^2) for Pseudo CN is:');
c = mean(const)


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

T = 1; %We set time a little higher to see if we can get better approximate order
a = 0.25;

trange = T./(5:2:11); %set of discretization intervals in time, i.e. the values of time(i)-time(i-1)
errorE = zeros(length(trange),1);
dx = sqrt(4*trange(1)/(2-a*trange(1))); %The minumum amplitude dx in order to keep stability

x = (-L:dx:L)';

count = 1;

%Explicit scheme
for dt = trange
    
    time = (0:dt:T);
    k = time(end); %Because in this way we are sure to compare the solutions at the same time
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*k);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*k);
    uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T
    
    N = length(x);
    r = dt/(2*dx^2);
        
    %Matrix for the explicit scheme
    d1 = (1-4*r)*ones(N,1);
    d2 = 2*r*ones(N-1,1);
    EM = diag(d1) + diag(d2,1) + diag(d2,-1);
    EM(1,1:2) = [1,0];
    EM(end,end-1:end) = [0,1];

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

const = errorE'./trange; %It gives the constant in time such that ||uNumerical - uExact||<= Const * Delta T.
disp('So the average constant for the O(Delta t) for Explicit is');
c = mean(const)

trange = T./trange; %Number of discretization points in the time domain

ord1 = (trange/trange(1)).^(-1)*errorE(1);
ord2 = (trange/trange(1)).^(-2)*errorE(1);
ord3 = (trange/trange(1)).^(-3)*errorE(1);
 
label = 'Tested number of points in time for dx='+string(dx);

figure;
loglog(trange,errorE,'r*',trange,ord1,'b-',trange,ord2,'k-',trange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in time measured at t=1 for Explicit scheme');
xlabel(label);
ylabel('Maximum error');

%Pseudo Crank-Nicolson
T = 30; %we fix a big time T so that we can have larger dt than dx^2 and the time error emerges correctly
dx = 0.2; %Small dx for the reason of beeing as precise as we can in space to let the error be almost fully depenedent on time
x = (-L:dx:L)';
N = length(x);

trange = T./[100,150,200,250];
errorC = zeros(length(trange),1);
count = 1;
for dt = trange
    time = (0:dt:T);
    r = dt/(2*dx^2);
    
    a1 = (1+2*r)*ones(N,1);
    a2 = -r*ones(N-1,1);
    A = diag(a1,0)+diag(a2,1)+diag(a2,-1);
    
    b1 = (1-2*r)*ones(N,1);
    b2 = r*ones(N-1,1);
    B = diag(b1,0) + diag(b2,1) + diag(b2,-1);
    
    m = 10^15;
    A(1,1) = m;
    A(end,end) = m;
    R = chol(A);
    
    k = time(end); %To be sure to compare the solution at the same time
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*k);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*k);
    uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T
    u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);
    u = u0;
        
    
    
    for t = time(2:end)
        b = B*u+dt*u.*(1-u).*(u-a);
        b(1) = 0;
        b(end) = m;
        v = R'\b;
        u = R\v;
    end
    errorC(count) = norm(u-uexact,inf);
    count = count + 1;
end

value = trange*a; %Even here we have all the values less than 2 so no problem with stability of Pseudo Crank Nicolson

const = errorC'./trange; %It gives the constant in time such that ||uNumerical - uExact||<= Const * Delta T.
disp('So the average constant for the O(Delta t) for Pseudo CN is:');
c = mean(const)
trange = T./trange;

ord1 = (trange/trange(1)).^(-1)*errorC(1);
ord2 = (trange/trange(1)).^(-2)*errorC(1);
ord3 = (trange/trange(1)).^(-3)*errorC(1);

label = 'Tested number of points in time for dx='+string(dx);

figure;
loglog(trange,errorC,'r*',trange,ord1,'b-',trange,ord2,'k-',trange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in time measured at t=30 for ps. Crank-Nicolson');
xlabel(label);
ylabel('Maximum error');

%So up to now we have compared the schemes in terms of accuracy which are the same. So
%since in terms of stability the pseuodo CN is less restrictive we should prefer it. 
%This is why we've chosen this scheme to do some tests varying the initial
%conditions in the other scripts.

% Just for a sake of completeness here it is the implementation of a numerical scheme
% which is second order both in space and time, it is not a finite differences scheme
% but after a finite difference discretization in space (of second order) we solve the
% pseudo discretized problem with a second order time integration, for instance a Runge Kutta 2.
% So what is actually implemented is METHOD OF LINES. This scheme has a
% bounded stability region too, so it is not the best possible choice but
% with these parameteres in space and time it shows the right convergence
% order.

a = 0.25;
T = 0.1;   
trange = 2.^(6:9);
N = 6002;
L = 100;
x = linspace(-L,L,N)';
dx = x(2)-x(1);

%Reference exact solution:
g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T

%Boundary conditions for this test are u(xL,t) = 0 and u(xR,t) = 1. 
    
% Initial condition for this test
u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);

d1 = -2*ones(N,1)/dx^2;
d2 = ones(N-1,1)/dx^2;
D2 = diag(d1,0) + diag(d2,-1) + diag(d2,1); %To define the stencil of centered finite differences
f = @(u) D2*u+u.*(1-u).*(u-a).*[0;ones(N-2,1);0];

error = zeros(length(trange),1);
count = 1;

for M = trange
    time = linspace(0,T,M)';
    u = u0;
    dt = time(2)-time(1);
    for t = time(2:end)
        ubar = u + dt*f(u); %Correction term with Euler scheme
        u = u + dt/2*(f(u)+f(ubar)); %Heun scheme
        u(1) = 0; %Boundary conditions
        u(end) = 1;

    end
    error(count) = norm(u-uexact,inf);
    count = count + 1;
end

ord1 = (trange/trange(1)).^(-1)*error(1);
ord2 = (trange/trange(1)).^(-2)*error(1);
ord3 = (trange/trange(1)).^(-3)*error(1);

figure;
loglog(trange,error,'r*',trange,ord1,'b-',trange,ord2,'k-',trange,ord3,'g-');
legend('error','order 1','order 2','order 3');
title('Convergence rate in time measured at t=0.1 for improved Euler method');
xlabel('Tested dt for dx=1e-4');
ylabel('Maximum error');
