clear all
close all
clc

%% Analysis of the travelling wave behaviour if time goes on

%Let's define a new setting, with Dirichlet BCs both equal to b, a number
%in [0,1] and see, varying the space interval, how the solution behaves as
%time goes on...in certain situations it will eventually converge to the
%stable stationary solutions u=0, u=a, u=1 and we can even show that this
%kind of convergence is faster for the Explicit method than for the Pseudo
%Crank Nicolson scheme

a = 0.25;
T = 50;   
M = 1001;
N = 301;
time = linspace(0,T,M);   
dt = time(2)-time(1);
H = 0.9; 

prompt = 'Which boundary condition u(-L,t)=u(L,t)=b do you want between b=1, b=0 or b=a?\n';
val = input(prompt);
switch val
    case 1
        b=1;
    case 0
        b=0;
    case a
        b=a;
    otherwise
        disp('Not an allowed input, let''s proceed with u(-L,t)=u(L,t)=a');
        b=a;
end
        

%Let's test the analytical results reported in some paper

criticalL = 0;
if b==0
    CriticalL = pi/(1-a);
    else if b==a
        CriticalL = pi;
        else if b==1
                CriticalL = pi/a;
            end
        end
end

%% Plot the surfaces checking the boundedness and asymptotic convergence
%% with Pseudo CN

if CriticalL == 0
    disp('Error: The test needs to be done with BCs wich are compatible with steady state solutions');
    else
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

        u1 = zeros(N,M);
        u2 = zeros(N,M);
        u1(:,1) = u0(x1,L1);
        u2(:,1) = u0(x2,L2);
        
        B1 = (speye(N)+CN1);
        B2 = (speye(N)+CN2);
        m = 10^15;
        
        %Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1
        B1(1,1) = m;
        B2(1,1) = m;
        B1(end,end) = m;
        B2(end,end) = m;
        R1 = chol(B1);
        R2 = chol(B2);
        
        i = 1;
        for t = time(2:end)
            b1 = (eye(N)-CN1)*u1(:,i)+dt*u1(:,i).*(1-u1(:,i)).*(u1(:,i)-a);
            b1(1) = m*b;
            b1(end) = m*b;
        
            v1 = R1'\b1;
            u1(:,i+1) = R1\v1;
            
            b2 = (eye(N)-CN2)*u2(:,i)+dt*u2(:,i).*(1-u2(:,i)).*(u2(:,i)-a);
            b2(1) = m*b;
            b2(end) = m*b;
        
            v2 = R2'\b2;
            u2(:,i+1) = R2\v2;
            
            i = i + 1;
        end
        
        figure;
        mesh(x1,time,u1')
        title('Solution in case L>CriticalL with a=0.25, b=1')
        figure;
        mesh(x2,time,u2')
        title('Solution in case L<CriticalL with a=0.25, b=1')
        
end

%% Check the "speed of convergence" of the method Pseudo CN against Explicit method

tolerance = 10.^(-8:1:-2);
L = CriticalL - 1;
x = linspace(-L,L,N)';
dx = x(2)-x(1);

s = 0.5; %So for stability we set dt = 2*s*dx^2/(4+a*dx^2)
% dt = 2*c*dx^2/(4+a*dx^2);
dt = 2*s*dx^2/(4+a*dx^2);

r = dt/(2*dx^2);
CN = 2*r*diag(ones(N,1))-r*diag(ones(N-1,1),1)-r*diag(ones(N-1,1),-1);

u0 = (b-H)*(x/L).^2+H; %So that it satisfies the BCs
u = u0;

%Matrix for the explicit scheme
d1 = (1-4*r)*ones(N,1);
d2 = 2*r*ones(N-1,1);
EM = diag(d1) + diag(d2,1) + diag(d2,-1);
EM(1,1:2) = [1,0];
EM(end,end-1:end) = [0,1];
runtimeEE = zeros(length(tolerance),1);
numberEE = zeros(length(tolerance),1);
count = 1;

for tol = tolerance
    tic;
    u = u0;
    err = 1;
    j=1;
    while err>tol
        u = EM*u + dt*u.*(1-u).*(u-a);
        j = j + 1;
        err = norm(u-b*ones(N,1),inf);
    end
    runtimeEE(count) = toc;
    numberEE(count) = j-1;
    count = count + 1;
end
figure;
semilogx(tolerance,numberEE,'r-*');
title('Number of iterations to converge to the steady solution');
legend('Number of iterations for explicit');
xlabel('Desired tolerances');
ylabel('Number of required iterations');


figure;
semilogx(tolerance,runtimeEE,'r-*');
title('Runtime required to reach convergence up to a given tolerance');
legend('Time required for the explicit');
xlabel('Desired tolerances');
ylabel('Time required for the execution');


%Same analysis for pseudo CN.
M = 1001;
N = 301;
% time = linspace(0,T,M);   
% dt = time(2)-time(1);

c=10^15;

B = (speye(N)+CN);

%Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1
B(1,1) = c;
B(end,end) = c;
R = chol(B);

runtimeCN = zeros(length(tolerance),1);
numberCN = zeros(length(tolerance),1);
count = 1;

for tol = tolerance
    tic;
    u = u0;
    err = 1;
    i=1;
    while err>tol
        f = (eye(N)-CN)*u+dt*u.*(1-u).*(u-a);
        f(1) = c*b;
        f(end) = c*b;

        u = R'\f;
        u = R\u;
        err = norm(u-b*ones(N,1),inf);
        i = i + 1;
    end
    numberCN(count) = i-1;
    runtimeCN(count) = toc;
    count = count + 1;
end

figure;
semilogx(tolerance,numberCN,'r-*');
title('Number of iterations to converge to the steady solution');
legend('Number of iterations for pseudo CN');
xlabel('Desired tolerances');
ylabel('Number of required iterations');


figure;
semilogx(tolerance,runtimeCN,'r-*');
title('Runtime required to reach convergence up to a given tolerance');
legend('Time required for pseudo CN');
xlabel('Desired tolerances');
ylabel('Time required for the execution');