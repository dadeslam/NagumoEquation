clear all
close all
clc

%% Analysis of the asymptotic behaviour towards constant stationary solutions

%Let's define a new setting, with Dirichlet BCs both equal to b, a number
%in {0,a,1} and see, varying the space interval, how the solution behaves as
%time goes on...in certain situations it will eventually converge to the
%stable stationary solutions u=0, u=a, u=1 and we can even show that this
%kind of convergence is faster for the Explicit method than for the Pseudo
%Crank Nicolson scheme

a = 0.25;
T = 50;   
M = 501;
N = 301;
time = linspace(0,T,M);   
dt = time(2)-time(1);
H = 0.5; %Parameter needed to define the initial condition

prompt = 'Which boundary condition u(-L,t)=u(L,t)=b do you want between b=1, b=0 or b=a?\n';
val = input(prompt);
switch val
    case 1
        b=1;
        CriticalL = pi/a;
    case 0
        b=0;
        CriticalL = pi/(1-a);
    case a
        b=a;
        CriticalL = pi;
    otherwise
        disp('Not an allowed input, let''s proceed with u(-L,t)=u(L,t)=a');
        b=a;
        CriticalL = pi;
end

%% Plot the surfaces checking the boundedness and asymptotic convergence
%% with Pseudo CN


%we just proceed in parallel with the 2 different domains to
%satisfy in one the bifurcation threshold and in other not
%satisfying it.
L1 = CriticalL + 1;
L2 = CriticalL - 1; %We should see a bifurcation in the 2 behaviours    
x1 = linspace(-L1,L1,N)';
x2 = linspace(-L2,L2,N)';
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);

r1 = dt/(2*dx1^2);
r2 = dt/(2*dx2^2);
CN1 = 2*r1*diag(ones(N,1))-r1*diag(ones(N-1,1),1)-r1*diag(ones(N-1,1),-1);
CN2 = 2*r2*diag(ones(N,1))-r2*diag(ones(N-1,1),1)-r2*diag(ones(N-1,1),-1);

u0 = @(y,L) (b-H)*(y/L).^2+H; %So that it satisfies the BCs

u1 = zeros(N,M); %I store all the solutions at any time
u2 = zeros(N,M);
u1(:,1) = u0(x1,L1);
u2(:,1) = u0(x2,L2);

B1 = (eye(N)+CN1);
B2 = (eye(N)+CN2);

m = 10^15;
%Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1 by keeping the SPD
%character of the matrices

B1(1,1) = m;
B2(1,1) = m;
B1(end,end) = m;
B2(end,end) = m;


R1=chol(B1);
R2=chol(B2);

i = 1;

for t = time(2:end)
    b1 = (eye(N)-CN1)*u1(:,i)+dt*u1(:,i).*(1-u1(:,i)).*(u1(:,i)-a);
    b1(1) = m*b;
    b1(end) = m*b;

    v = R1'\b1;
    u1(:,i+1) = R1\v;

    b2 = (eye(N)-CN2)*u2(:,i)+dt*u2(:,i).*(1-u2(:,i)).*(u2(:,i)-a);
    b2(1) = m*b;
    b2(end) = m*b;

    v = R2'\b2;
    u2(:,i+1) = R2\v;
    
    i = i + 1;
end

figure;
mesh(x1,time,u1')
title('Solution in case L>CriticalL with a=0.25, b='+string(b))
figure;
mesh(x2,time,u2')
title('Solution in case L<CriticalL with a=0.25, b='+string(b))
        


%% Check the "speed of convergence" of the method Pseudo CN against Explicit method with the case b=a=0.25, 
%% testing both H = 0.1 and H=0.3

b=a;
N = 5;
ic = [0.1,0.3];
CriticalL = pi;
tolerance = 10.^(-4:1:-1); %Tested tolerance
L = CriticalL - 1; %I want to get asymptotic convergence so I impose the sufficient condition to be satisfied
x = linspace(-L,L,N)';
dx = x(2)-x(1);

s = 0.99; %So for stability we set dt = 2*s*dx^2/(4+a*dx^2)
dt = 2*s*dx^2/(4+a*dx^2);

for H = ic
    u0 = (b-H)*(x/L).^2+H; %So that it satisfies the BCs
    u = u0;
    r = dt/(2*dx^2);

    %Matrix for the explicit scheme
    d1 = (1-4*r)*ones(N,1);
    d2 = 2*r*ones(N-1,1);
    EM = diag(d1) + diag(d2,1) + diag(d2,-1); %Matrix for the Explicit scheme
    EM(1,1:2) = [1,0];
    EM(end,end-1:end) = [0,1];
    numberEE = zeros(length(tolerance),1); %Here I will store the number of iterations required to reach those tolerances
    count = 1;

    for tol = tolerance
        u = u0;
        err = 1;
        j=1;
        while err>tol
            u = EM*u + dt*u.*(1-u).*(u-a);
            j = j + 1;
            err = norm(u-b,inf);
        end
        numberEE(count) = j-1; %Because I've increased j one time more at the end
        count = count + 1;
    end
    uE = u;

    figure;
    tit = 'Number of iterations to converge to the steady solution b=a and for H='+string(H);
    semilogx(tolerance,numberEE,'r-*');
    title(tit);
    legend('Number of iterations for explicit');
    xlabel('Desired tolerances');
    ylabel('Number of required iterations');

    %Same analysis for pseudo CN.
    r = dt/(2*dx^2);
    CN = 2*r*diag(ones(N,1))-r*diag(ones(N-1,1),1)-r*diag(ones(N-1,1),-1);

    B = (eye(N)+CN); %Matrix for the Pseudo CN method

    B(1,1:2) = [r,0];
    B(end,end-1:end) = [0,r];

    numberCN = zeros(length(tolerance),1);
    count = 1;
    for tol = tolerance
        u = u0;
        err = 1;
        i=1;
        while err>tol
            f = (eye(N)-CN)*u+dt*u.*(1-u).*(u-a);
            f(1) = r*b;
            f(end) = r*b;
            u = B\f;
            err = norm(u-b,inf);
            i = i + 1;
        end
        numberCN(count) = i-1;
        count = count + 1;
    end
    uCN = u;
    figure;
    tit = 'Number of iterations to converge to the steady solution for b=a and H='+string(H);
    semilogx(tolerance,numberCN,'r-*');
    title(tit);
    legend('Number of iterations for pseudo CN');
    xlabel('Desired tolerances');
    ylabel('Number of required iterations');
end