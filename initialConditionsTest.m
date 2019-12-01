%% Let's now analyze various initiali conditions, using the Pseudo Crank Nicolson method because we prefer not
%% having restrictions on the discretization steps to get stability, so that we can play
%% with those values as much as we want
clc;
clear all;
close all;

global x;
global time;
global a;

xL = -100;
xR = 100;
T = 30;
a = 0.25;

M = 1001;
N = 501;


time = linspace(0,T,M);   
x = linspace(xL,xR,N)';

dx = x(2)-x(1);
dt = time(2)-time(1);

r = dt/(2*dx^2);
a1 = 2*r*ones(N,1);
a2 = -r*ones(N-1,1);
A1 = diag(a1);
A2 = diag(a2,1);
CN = A1+A2+A2';
% Now we will work with Dirichlet homogeneus BCs

L = 30; %Parameter used to define step functions

%% Initial conditions with single connected non zero region 

prompt = 'In the initial condition (x>-L).*(x<L)*c, do you want c<a, c=a or c>a?\n Type 1 for the first choice\n Type 2 for the second one \n Type 3 for the third\n';
val = input(prompt);
switch val
    case 1
        u0 = (x>-L).*(x<L)*a/2; %Impulse less than a 
    case 2
        u0 = (x>-L).*(x<L)*a; %Impulse equal to a 
    case 3
        u0 = (x>-L).*(x<L)*3*a/2; %Impulse greater than a
    otherwise
        disp('Not an admissible choice, let''s assume you have chosen the third one')
        u0 = (x>-L).*(x<L)*3*a/2; %Impulse greater than a
end

u = solve(u0,CN);
figure;
mesh(x,time,u')
title('Solution with a=0.25')
xlabel('Space domain');
ylabel('Time domain');
zlabel('Value of the solution u=u(x,t)')


%% Initial conditions with 2 connected and disjoint subintervals of [xL,xR]different from 0

prompt = 'In the initial condition (x>-2L).*(x<-L)*c+(x>L).*(x<2L)*d, do you want c,d<a, c>a & d<a or c,d>a?\n Type 1 for the first choice\n Type 2 for the second one \n Type 3 for the third\n';
val = input(prompt);
switch val
    case 1
        u0 = (x>-2*L).*(x<-L)*2*a/3 + (x>L).*(x<2*L)*2*a/3; %Sum of the 2 impulses
        %greater than a but both alone less than a 
    case 2
        u0 = (x>-2*L).*(x<-L)*5*a/3 + (x>L).*(x<2*L)*a/3; %Sum of the 2 impulses greater than a and one of them greater than a alone
    case 3
        u0 = (x>-2*L).*(x<-dx)*5*a/3 + (x>dx).*(x<2*L)*4*a/3; %Both greater than a and quite near, the solution smooths out as time goes on
        %Important thing: It is not important how big is the sum of the local
        %maxima in order to see if the solution will attain value 1 or 0 as time
        %goes on...the important thing is that locally exists some point where u0
        %is greater sctrictly than a
    otherwise
        disp('Not an admissible choice, let''s assume you have chosen the third one')
        u0 = (x>-2*L).*(x<-dx)*5*a/3 + (x>dx).*(x<2*L)*4*a/3;
end
u = solve(u0,CN);
figure;
mesh(x,time,u')
title('Solution with a=0.25')
xlabel('Space domain');
ylabel('Time domain');
zlabel('Value of the solution u=u(x,t)')

%% Smooth initial condition with just a little number of discretization points with values over a

prompt = 'In the initial condition c*e^(-x.^2./50), do you want c=a,c a little greater than a or c "a lot" greater than a?\n Type 1 for the first choice\n Type 2 for the second one \n Type 3 for the third\n';
val = input(prompt);
switch val
    case 1
        u0 = exp(-x.^2./50)*a;    
    case 2
        u0 = exp(-x.^2./50)*(12/10*a); %even if on 21 discretization points u0>a, the solution will converge to 0
        % %so there is a sort of minimum "power" to give to the system in order to
        % %converge to 1 as time goes on 
        disp('Here the number of discretization points where u0(x)>a is:\n')
        n_points_u0_overA = sum(u0>a)
    case 3
       u0 = exp(-x.^2./50)*(14/10*a); %here with 29 point of u0 over a I converge to 1 as time goes on
       disp('Here the number of discretization points where u0(x)>a is:\n')
       n_points_u0_overA = sum(u0>a) %---> 9
    otherwise
        disp('Not an admissible choice, let''s assume you have chosen the third one')
        u0 = exp(-x.^2./50)*(14/10*a); %here with 29 point of u0 over a I converge to 1 as time goes on
end
u = solve(u0,CN);
figure;
mesh(x,time,u')
title('Solution with a=0.25')
xlabel('Space domain');
ylabel('Time domain');
zlabel('Value of the solution u=u(x,t)')


function u = solve(u0,CN)
    % The scheme writes, in matrix form, as follows:
    % (I+A)u_n+1 = (I-A)*u+dt*u.*(1-u).*(u-a)
    %%Plot of the surface representing the solution to the PDE
    
    global x;
    global time;
    global a;
    
    N = length(x);
    M = length(time);
    dt = time(2)-time(1);
    u = zeros(N,M);
    u(:,1) = u0;

    B = (speye(N)+CN);
    m = 10^15;
    
    %Setting for the right BCs: u(xL,t) = 0, u(xR,t) = 1
    B(1,1) = m;
    B(end,end) = m;
    R = chol(B);
    
    i = 1;
    for t = time(2:end)
        b = (eye(N)-CN)*u(:,i)+dt*u(:,i).*(1-u(:,i)).*(u(:,i)-a);
        b(1) = 0;
        b(end) = 0;

        v = R'\b;
        u(:,i+1) = R\v;
        i = i + 1;
    end
end