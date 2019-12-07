clear all
close all
clc

%Solution of Nagumo equation with a piecewise linear finite elements'
%aproach

%% Guide on the structure of the script:
%% The finite elements method is implemented in the function FEM which gets as
%% inputs the initial condition. Then all the parameters are set inside the function,
%% while there are 2 global objects which are x and a which are accessible to that function too

%% Then in the main there is the sequence of the 3 switch cases where 3 families of 
%% initial conditions are tested like in the finite differences case.

global x %I set these global so that I can access to them in the function
global a

a = 0.25;
L = 100;

N = 501; %Tested number of points to discretize the space domain 
x = linspace(-L,L,N)';
dx = x(2)-x(1);


%% Initial conditions with single connected non zero region 
Lbar = 30;
prompt = 'In the initial condition (x>-L).*(x<L)*c, do you want c<a, c=a or c>a?\n Type 1 for the first choice\n Type 2 for the second one \n Type 3 for the third\n';
val = input(prompt);
switch val
    case 1
        u0 = (x>-Lbar).*(x<Lbar)*a/2; %Impulse less than a 
    case 2
        u0 = (x>-Lbar).*(x<Lbar)*a; %Impulse equal to a 
    case 3
        u0 = (x>-Lbar).*(x<Lbar)*3*a/2; %Impulse greater than a
    otherwise
        disp('Not an admissible choice, let''s assume you have chosen the third one')
        u0 = (x>-Lbar).*(x<Lbar)*3*a/2; %Impulse greater than a
end

figure;
[u,time] = FEM(u0);
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
        u0 = (x>-2*Lbar).*(x<-Lbar)*2*a/3 + (x>Lbar).*(x<2*Lbar)*2*a/3; %Sum of the 2 impulses
        %greater than a but both alone less than a 
    case 2
        u0 = (x>-2*Lbar).*(x<-Lbar)*5*a/3 + (x>Lbar).*(x<2*Lbar)*a/3; %Sum of the 2 impulses greater than a and one of them greater than a alone
    case 3
        u0 = (x>-2*Lbar).*(x<-dx)*5*a/3 + (x>dx).*(x<2*Lbar)*4*a/3; %Both greater than a and quite near, the solution smooths out as time goes on
        %Important thing: It is not important how big is the sum of the local
        %maxima in order to see if the solution will attain value 1 or 0 as time
        %goes on...the important thing is that locally exists some point where u0
        %is greater sctrictly than a
    otherwise
        disp('Not an admissible choice, let''s assume you have chosen the third one')
        u0 = (x>-2*Lbar).*(x<-dx)*5*a/3 + (x>dx).*(x<2*Lbar)*4*a/3;
end
figure;
[u,time] = FEM(u0);
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
        disp('Here the number of discretization points where u0(x)>a is:')
        n_points_u0_overA = sum(u0>a)
    case 3
       u0 = exp(-x.^2./50)*(14/10*a); %here with 29 point of u0 over a I converge to 1 as time goes on
       disp('Here the number of discretization points where u0(x)>a is:')
       n_points_u0_overA = sum(u0>a) %---> 9
    otherwise
        disp('Not an admissible choice, let''s assume you have chosen the third one')
        u0 = exp(-x.^2./50)*(14/10*a); %here with 29 point of u0 over a I converge to 1 as time goes on
end
figure;
[u,time] = FEM(u0);
mesh(x,time,u')
title('Solution with a=0.25')
xlabel('Space domain');
ylabel('Time domain');
zlabel('Value of the solution u=u(x,t)')



function [u,time] = FEM(u0)
    
    global x;
    global a;

    T = 30;
    N = length(x);
    M = 400; %Number of time steps
    time = linspace(0,T,M); %Time discretized domain
    dt = time(2)-time(1); %Time interval

    h = diff(x); %I compute all the interval amplitudes x(i+1)-x(i) in this way

    %Assembling the stiffness matrix:
    d1 = [1/h(1);1./h(1:end-1) + 1./h(2:end);1/h(end)];
    d2 = - 1./h(1:end);
    A = diag(d1,0) + diag(d2,1) + diag(d2,-1); 

    %Assembling the mass matrix
    m1 = [h(1)/3;1/3*(h(1:end-1)+h(2:end));h(end)/3];
    m2 = 1/6*h; %Upper diagonal
    Mass = diag(m1,0) + diag(m2,1) + diag(m2,-1);

    index = 2:N-1;  

    Left = (Mass+dt/2*A); %Matrices involved in the loop
    Right = (Mass-dt/2*A);

    c = 10^15;

    Left(1,1) = c; 
    Left(end,end)= c; %penalization method to keep simmetry of left matrix and use cholesky factorization
    R = chol(Left);
    l=0;
    r=0;

    u = zeros(N,M);
    u(:,1) = u0;
    count = 1;
    for t = 1:M-1
        F = zeros(N,1);
        %F(index) = RHS(u,index); %This calls the function computing the integral using Cavalieri-Simpson scheme
        F(index) = gauss(u(:,count),index); %This calls the function computing the integral using Gauss quadrature
        v = R'\(Right*u(:,count)+dt*F+[c*l;zeros(N-2,1);c*r]);
        u(:,count+1) = R\v;
        count = count + 1;
    end
end

%Function computing the integral of f(u)*phi_i(x) with Cavalieri-Simpson
%rule
function F = RHS(u,i)
    global x
    global a
    h=diff(x);
    
    F = h(i-1)/6*4*1/2.*((u(i)/2+u(i-1)/2).*(1-u(i-1)/2-u(i)/2).*(u(i)/2+u(i-1)/2-a))+...
        h(i-1)/6 .* u(i) .* (1-u(i)) .* (u(i)-a) + ...
        h(i)/6 * 4 *1/2 .* ((u(i)/2+u(i+1)/2).*(1-u(i+1)/2-u(i)/2).*(u(i)/2+u(i+1)/2-a))+...
        h(i)/6 .* u(i).* (1-u(i)) .* (u(i)-a);
    
end


%Implementation of Gaussian quadrature with 3 points, precise for 2*3-1=5>4
%so for our case it is:
function F = gauss(u,i)
    global x
    global a

    phi0 = @(y) (y-x(i-1))./(x(i)-x(i-1)).*(y>x(i-1)).*(y<x(i)) + (x(i+1)-y)./(x(i+1)-x(i)).*(y>=x(i)).*(y<x(i+1));%phi_i
    phil = @(y) (x(i)-y)./(x(i)-x(i-1)).*(y>=x(i-1)).*(y<x(i)); %phi_i+1 restricted to (x_i-1,x_i)
    phir = @(y) (y-x(i))./(x(i+1)-x(i)).*(y>x(i)).*(y<=x(i+1)); %phi_i-1 restricted to (x_i,x_i+1)
    vl = @(y) u(i-1).*phil(y) + u(i).*phi0(y); %Is the u in the interval [x_i-1,x_i]
    vr = @(y) u(i).*phi0(y) + u(i+1).*phir(y); %Is the u in the interval [x_i,x_i+1]
    mapr = @(y) 0.5*(1-y).*x(i)+0.5*(1+y).*x(i+1); %Change of variable from x(i),x(i+1) to -1,1.
    mapl = @(y) 0.5*(1-y).*x(i-1)+0.5*(1+y).*x(i); %Change of variable from x(i-1),x(i) to -1,1.
    fl = @(y) phil(mapl(y)).*vl(mapl(y)).*(1-vl(mapl(y))).*(vl(mapl(y))-a);
    fr = @(y) phir(mapr(y)).*vr(mapr(y)).*(1-vr(mapr(y))).*(vr(mapr(y))-a);

    %fl = @(y) vl(mapl(y)); %These are just two lines which if uncommented
                            %allow to integrate the function discretized by the vector u, just to
                            %test that the quadrature rule works.
    %fr = @(y) vr(mapr(y));

    Fl = (5/9*(fl(-sqrt(3/5))+fl(sqrt(3/5)))+8/9*fl(0)).*(x(i)-x(i-1))/2;
    Fr = (5/9*(fr(-sqrt(3/5))+fr(sqrt(3/5)))+8/9*fr(0)).*(x(i+1)-x(i))/2;
    F = Fl + Fr;
end