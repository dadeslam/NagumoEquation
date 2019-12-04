clear all
close all

global x
global a
a = 0.25;

%% Speed of convergence to the asymptotic stationary solution


T = 50;  %Final time, never reached but just to give a time step to advance with in the execution  
M = 1001; %Number of points in the interval [0,T]
N = 501; %Number of space points
time = linspace(0,T,M);   
dt = time(2)-time(1);
H = 0.9;  %Parameter used to define the initial condition
b = a; %Value taken by the solution at both the space boundaries : +L and -L

if b==0 %Defining the bifurcation lengths according to the analytical study in the report
    CriticalL = pi/(1-a);
    else if b==a
        CriticalL = pi;
        else if b==1
                CriticalL = pi/a;
            end
        end
end

%% Check the speed of convergence of the method to the stationary
% constant solution

tolerance = 10.^(-1:-1:-4); %Vector of fixed desired tolerances to reach
L = CriticalL - 1; %I want that the sufficient condition is satisfied
x = linspace(-L,L,N)';
h = diff(x); %I find the discretization intervals x(i)-x(i-1)

u0 = (b-H)*(x/L).^2+H; %So that it satisfies the BCs

index = 2:N-1;  %Indices of hat functions which define the basis of my space of piecewise linear functions

%Assembling the stiffness matrix:
d1 = [1/h(1);1./h(1:end-1) + 1./h(2:end);1/h(end)];
d2 = - 1./h(1:end);
A = diag(d1,0) + diag(d2,1) + diag(d2,-1);

%Assembling the mass matrix
m1 = [h(1)/3;1/3*(h(1:end-1)+h(2:end));h(end)/3];
m2 = 1/6*h; 
Mass = diag(m1,0) + diag(m2,1) + diag(m2,-1);


Left = (Mass+dt/2*A); %Matrices involved in the loop to solve the problem
Right = (Mass-dt/2*A);

c = 10^15; %Penalization parameter

Left(1,1) = c;
Left(end,end)= c; %penalization method to keep simmetry of left matrix and use cholesky factorization
R = chol(Left);

number = zeros(length(tolerance),1);
runtime = zeros(length(tolerance),1);
count = 1;

for tol = tolerance
    tic; %Start measuring time
    u = u0;
    err = 1; %Set it higher than any tolerance to get inside the while loop
    i=1;
    while err>tol
        F = zeros(N,1);
        F(index) = RHS(u,index); %Computing the right hand side with Cavalieri-Simpson
        u = R'\(Right*u+dt*F+[c*b;zeros(N-2,1);c*b]);
        u = R\u;
        err = norm(u-b*ones(N,1),inf); %Computing the infinity norm of the error
        i = i + 1;
    end
    number(count) = i-1;
    runtime(count) = toc; %Stop measuring time and storing the result
    count = count + 1;
end

figure;
semilogx(tolerance,number,'r-*');
title('Number of iterations required to get the desired convergence');
xlabel('Desired tolerances');
ylabel('Number of required iterations');


figure;
semilogx(tolerance,runtime,'r-*');
title('Runtime required to reach convergence up to a given tolerance');
xlabel('Desired tolerances');
ylabel('Time required for the execution');

%Function computing the integral of f(u)*phi_i(x) with Cavalieri-Simpson
function F = RHS(u,i)
    global x
    global a
    h=diff(x);
    
    F = h(i-1)/6*4*1/2.*((u(i)/2+u(i-1)/2).*(1-u(i-1)/2-u(i)/2).*(u(i)/2+u(i-1)/2-a))+...
        h(i-1)/6 .* u(i) .* (1-u(i)) .* (u(i)-a) + ...
        h(i)/6 * 4 *1/2 .* ((u(i)/2+u(i+1)/2).*(1-u(i+1)/2-u(i)/2).*(u(i)/2+u(i+1)/2-a))+...
        h(i)/6 .* u(i).* (1-u(i)) .* (u(i)-a);
    
end