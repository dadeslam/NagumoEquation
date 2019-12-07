clear all
close all


%Solution of Nagumo equation with a piecewise linear finite elements'
%aproach

L = 100; %space delmiter
T = 0.1; %time delimiter

global x %I set these global so that I can access to them in the function
global a

M = 400; %Number of time steps
a = 0.25; %Threshold parameter
time = linspace(0,T,M); %Time discretized domain
dt = time(2)-time(1); %Time interval

xrange = (400:400:2000); %Tested number of points to discretize the space domain
%xrange = xrange(2); %used to get the plot in the report --->uncomment this
%line to get the plot in the report, and uncomment lines 99-105 below.
error = zeros(length(xrange),1); %will be used to store the infinity norm of the errors
E = zeros(length(xrange),1); %will be used to store the H1 semi-norms of the errors
count = 1;
for N = xrange
    
    N1 = floor(N/6); %number of space steps in between -100 and -45, and even in 45 e 100
    N2 = floor(4*N/6); %number of space steps in the interval -45 e 45 (just a way to refine a non uniform mesh)
    x1 = linspace(-L,-45,N1);
    x2 = linspace(-45,45,N2);
    x3 = linspace(45,L,N1);
    x = [x1,x2(2:end-1),x3]';
    N = length(x); %Since the refining strategy is not perfect at the end we will get a number of points which is slightly different than original N.
    xrange(count) = N; %So I update the value in the original vector of points.
    
    h = diff(x); %I compute all the interval amplitudes x(i+1)-x(i) in this way
    
    %Assembling the stiffness matrix:
    d1 = [1/h(1);1./h(1:end-1) + 1./h(2:end);1/h(end)];
    d2 = - 1./h(1:end);
    A = diag(d1,0) + diag(d2,1) + diag(d2,-1);

    %Assembling the mass matrix
    m1 = [h(1)/3;1/3*(h(1:end-1)+h(2:end));h(end)/3];
    m2 = 1/6*h; %Upper diagonal
    Mass = diag(m1,0) + diag(m2,1) + diag(m2,-1);

    %Reference exact solution:
    g1 = 0.5*(sqrt(2)*x+(1-2*a)*T);
    g2 = 0.5*(sqrt(2)*a*x+a*(a-2)*T);
    uexact = (exp(g1)+a*exp(g2))./(exp(g1)+exp(g2)+1); %Exact solution at time T
    u0 = (exp(sqrt(2)*x/2)+a*exp(sqrt(2)*a*x/2))./(exp(sqrt(2)*x/2)+exp(sqrt(2)*a*x/2)+1);

    l=0; %Left and right Dirichlet boundary conditions
    r=1;

    index = 2:N-1;  

    Left = (Mass+dt/2*A); %Matrices involved in the loop
    Right = (Mass-dt/2*A);

    c = 10^15;

    Left(1,1) = c; 
    Left(end,end)= c; %penalization method to keep simmetry of left matrix and use cholesky factorization
    R = chol(Left);

    u = u0;

    for t = 1:M-1
        F = zeros(N,1);
        F(index) = RHS(u,index); %This calls the function computing the integral using Cavalieri-Simpson scheme
        %F(index) = gauss(u,index); %This calls the function computing the integral using Gauss quadrature
        u = R'\(Right*u+dt*F+[c*l;zeros(N-2,1);c*r]);
        u = R\u;
    end
    E(count) = energy(u-uexact,A); %Compute the H1 semi-norm of the error
    error(count) = norm(u-uexact,inf); %Compute the infinity norm of the error
    count = count + 1;
end

disp('H1 seminorm of the difference decreasing as the space discretization is refined, t=0.1')
E
disp('Infinity norm of the difference decreasing as the space discretization is refined,t=0.1')
error

figure;
semilogy(xrange,E,'r-*');
title('Decay in the H1 semi-norm of the error');
xlabel('Number of discretization points');
ylabel('Error in the H1 semi-norm');

figure;
semilogy(xrange,error,'r-*');
title('Decay in the infinity norm of the error');
xlabel('Number of discretization points');
ylabel('Error in the infinity norm');

% tit = 'Comparison of the solutions';
% %case with 800 space points.
% plot(x,u,'ro',x,uexact,'k','Markersize',4);
% title(tit);
% legend('Numerical solution','Exact solution');
% xlabel('Space coordinates');
% ylabel('Solution u(x,T)');

%Function computing the H1 semi-norm of the vector v, just by using matrix
%A which is the stiffness matrix.
function e = energy(v,A)
    global x;
    e = 0;
    for i = 1:length(x)
        for j = 1:length(x)
            e = e + v(i)*v(j)*A(i,j);
        end
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