%% MHA 021 FEM
% Assignment 1
% Group 33
% Nils Helgesson, Gabriel Wendel
%% Task 1
close all;
clear all;
clc;

%% 1D heat
% Indata
A=1;                    % Area [m^2]   
alpha=1;                % Convection heat transfer coefficient (same for inner and outer)
T_in=20;                % Temperature inside [degrees C]
T_out=7.3;              % Temperature outside [degrees C]

% Layer thickness [m]
h1=50e-3;
h2=250e-3;
h3=2e-3;
h4=25e-3;
h5=12e-3;

% Array of layers
h=[h5 h4 h3 h2 h1];

% Heat transfer coefficient for all materials [W/m^2/K]
k1=15;
k2=0.04;
k3=5;
k4=k2;
k5=0.25;

% Array of heat transfer coefficients
k=[k5 k4 k3 k2 k1];

% Define x-coordinates
L1=h(1)
L2=sum(h(1:2))
L3=sum(h(1:3))
L4=sum(h(1:4))
L=sum(h(1:5))

coords=[0 L1 L2 L3 L4 L];
nnodes=length(coords);                  % Number of nodes
nel=nnodes-1;                           % Number of elements
dofs=nnodes;                            % Number of Degrees of freedom

Edof=[(1:nel); (1:nel); (2:nnodes)]';   % Degrees of freedom for elements

K=zeros(dofs);                          % Initialize stiffness matrix      

% Initialize shape function matrices
N_L=zeros(1,nnodes);
N_L(end)=1;                             % Last shape function, N(L) = 1
N_0=zeros(1,nnodes);
N_0(1)=1;                               % First shape function, N(0) = 1


for i=1:nel
    syms x
    x_1=coords(i);
    x_2=coords(i+1);
    
    L_e=x_2-x_1;
    % N_e = [N_1^e N_2^e]
    N_e=[-1/L_e*(x-x_2), 1/L_e*(x-x_1)];
    B_e=diff(N_e,x);
    % Stiffness matrix for element
    Ke=int(B_e'*k(i)*A*B_e,x,x_1,x_2);
    % Assemble using CALFEM function "assem"
    K=assem(Edof(i,:),K,Ke);
end

% stiffness matrix (convective)
Kc=N_L'*A*alpha*N_L+N_0'*A*alpha*N_0;

K=K+Kc

% load vector (convective)
fc=N_L'*A*alpha*T_out+N_0'*A*alpha*T_in

% solve to obtain vector a
a=solveq(K,fc)

% plot result
figure(1)
plot(coords,a,'o-','LineWidth',1.9)
xlabel('Distance trough roof x[m]')
ylabel('Temperature [{\circ} C]')
title('Temperature distribution')
hold on

% Plot analytical solution of temperature distribution through wall with
% constant heat conductivity
k_mean = 1/(sum(h ./ k));
x = linspace(0, L, 100); % coordinate through the wall
T = (T_in*alpha*k_mean + T_out*alpha*k_mean + L*T_in*alpha*alpha -...
T_in*alpha*alpha .* x + ...
T_out*alpha*alpha .* x)/(alpha*k_mean + alpha*k_mean + L*alpha*alpha);
plot(x, T,'LineWidth',1.9)
legend('FEM solution','Analytical Solution')



% Heat loss per square meter
% q=zeros(1,nel);
% figure(2)
% hold on
% for i=1:nel
%     x_1=coords(i);
%     x_2=coords(i+1);
%     L_e=x_2-x_1;
%     T_i=a(i)
%     T_j=a(i+1)
%     q(i)=-k(i)*(T_j-T_i)/L_e
%     plot([x_1 x_2],[q(i) q(i)],'o-','LineWidth',1.9)
% end

dt_dx = zeros(1,nel);
figure(2)
hold on

% FE-approximation
for i=1:nel
    syms x
    x_1=coords(i);
    x_2=coords(i+1);
    L_e=x_2-x_1;    % Element Length
    a_e = [a(i); a(i+1)];
    % N_e = [N_1^e N_2^e]
    N_e=[-1/L_e*(x-x_2), 1/L_e*(x-x_1)];
    B_e=diff(N_e,x);    % B = dN/dx
    dt_dx(i) = -k(i)*B_e*a_e;   % du/dt = Ba
    plot([x_1 x_2],[dt_dx(i), dt_dx(i)],'o-','LineWidth',1.9)
end

ylim([1.42 1.43])
xlabel('Distance trough roof x[m]')
ylabel('Heat loss per square meter q [{W/m^2}]')



%% Task 2
clear all
clc
%Variables to change
water_depth=100;        % Nordstream mean depth for an example [m]
nel=5;                  % Number of elements

%indata
rho=1000;               % Density [kg/m^3]
g=9.82;                 % Gravitational constant [m/s^2]
p_atm=101300;           % Atmospheric pressure [Pa]
a=0.425;                % Inner radius [m]
b=0.450;                % Outer radius [m]
p_i=10e6;               % Pressure inside pipe [Pa]
E=210e9;                % Youngs modulus [Pa]
nu=0.3;                 % Poisson's ratio


p_e=p_atm+rho*g*water_depth;    % Total pressure


nnodes=nel+1;                   % Number of nodes
dofs=nnodes;                    % Degrees of freedom
coords=linspace(a,b,nnodes);    % Coordinates
h=(b-a)/(nnodes-1);             % Length of elements

Edof=[(1:nel); (1:nel); (2:nnodes)]';   % Elements degrees of freedom

% Initialize stiffness matrix
K=zeros(dofs);

for i=1:nel
    r_i=coords(i);
    r_ip1=coords(i+1);
    syms r

    % Comparison, used the expression for element stiffness matrix (Task c))
    % N1e=(r_ip1-r)/h;
    % N2e=(r-r_i)/h;
    % B1e=diff(N1e,r);
    % B2e=diff(N2e,r);
    % f_11=N1e^2/r+2*N1e*B1e+r*B1e^2;
    % f_12=N1e*N2e/r+N1e*B2e+N2e*B1e+r*B1e*B2e;
    % f_21=f_12;
    % f_22=N2e^2/r+2*N2e*B2e+r*B2e^2;
    % Ke=E/(1-nu^2)*int([f_11,f_12;f_21,f_22],r,r_i,r_ip1);
    
    % Calculate components of element stiffness matrix according to the
    % given expression for K^e
    f_11=-2+(r_ip1/h)^2*log(r_ip1/r_i);
    f_12=-r_i*r_ip1/h^2*log(r_ip1/r_i);
    f_21=f_12;
    f_22=2+(r_i/h)^2*log(r_ip1/r_i);
    % Elements stiffness matrix
    Ke=E/(1-nu^2)*[f_11,f_12;f_21,f_22];
    % Assemble
    K=assem(Edof(i,:),K,Ke);
end

% Initialize shape function matrices
N_b=zeros(1,nnodes);
N_b(end)=1;
N_a=zeros(1,nnodes);
N_a(1)=1;

% Stiffness matrix
Kc=-(E/(1+nu))*(N_b'*N_b-N_a'*N_a);

% Load vector
fc=a*N_a'*p_i-b*N_b'*p_e;

K=K+Kc;

% Solve for deformation, u
u=solveq(K,fc);

% Plot results
figure(1)
plot(coords,u,'ro-','LineWidth',1.9)
hold on

% Plot and compare with analytical solution

% Derived expressions for A, B, A1 and A2, see report
A=(b^2/(b^2-a^2))*(p_i-p_e)-p_i;
B=(a^2*b^2/(b^2-a^2))*(p_i-p_e);
A1=A*(1-nu)/E;
A2=B*(1+nu)/E;

u_analytic=A1*r+A2/r;
figure(1)
fplot(r,u_analytic,[a b],'y--','LineWidth',1.9)
legend('u(r) FE solution','u(r) analytic solution')
xlabel('r [m]')
ylabel('u(r) [m]')

% Stresses
sigma_r_analytic=A-B/r^2;
sigma_theta_analytic=A+B/r^2;

% Compute stresses at the center of each element
% for i=1:nel
%     % Derivative of element deformation
%     du_dr=(u(i+1)-u(i))/h;
%     % deformation at center
%     u_center=(u(i+1)+u(i))/2;
%     r_center(i)=(coords(i)+coords(i+1))/2;
%     % Compute stresses according to Eq (2.3)
%     sigma_r_FE(i)=(E/(1-nu^2))*(du_dr+nu*u_center/r_center(i));
%     sigma_theta_FE(i)=(E/(1-nu^2))*(u_center/r_center(i)+nu*du_dr);
% end

% Be = 1/h*[-1 1];    % Be = dN/dr

%r_center_local = h/2;     % center of element for local coordinates

sigma_r_FE = zeros(1,nel);
sigma_theta_FE = zeros(1,nel);
r_center = zeros(1,nel);

N_r = [0.5, 0.5]; % N_r = [1-r_center(i)/h r_center(i)/h] where r_center = Le/2
Be = 1/h*[-1, 1]; % B = dN/dx


for i=1:nel
    r_center(i)=(coords(i)+coords(i+1))/2;
    % Compute stresses according to Eq (2.3)
    sigma_r_FE(i)=(E/(1-nu^2))*(Be*[u(i); u(i+1)]+nu*(N_r*[u(i); u(i+1)])/r_center(i));
    sigma_theta_FE(i)=(E/(1-nu^2))*(N_r*[u(i); u(i+1)]/r_center(i)+nu*(Be*[u(i); u(i+1)]));
end




% Plot stresses
figure(2)
plot(r_center,sigma_r_FE,'ro-','LineWidth',1.9)
hold on
fplot(r,sigma_r_analytic,[a b],'y--','LineWidth',1.9)
xlabel('r [m]')
ylabel('{\sigma_r(r) [Pa]}')
legend('{\sigma_r(r)} FE','{\sigma_r(r)} analytic','Location','southeast')
title('Comparison of radial stresses')

figure(3)
plot(r_center,sigma_theta_FE,'ro-','LineWidth',1.9)
hold on
fplot(r,sigma_theta_analytic,[a b],'y--','LineWidth',1.9)
xlabel('r [m]')
ylabel('{\sigma_\theta(r) [Pa]}')
legend('{\sigma_\theta(r)} FE','{\sigma_\theta(r)} analytic')
title('Comparison of tangential stresses')

% Find max value for deformation and stresses
u_max=max(u)
sigma_r_max=max(sigma_r_FE)
sigma_theta_max=max(sigma_theta_FE)



%% Convergence analysis
% We choose to compare max deformation and max stresses by running the code
% multiple times in order to find an appropriate number of elements.

% Number of elements
nel_conv=[1 2 4 8 16 32 64 128];
% Resulting max deformation
u_max=[3.1989 3.1993 3.1993 3.1994 3.1994 3.1994 3.1994 3.1994]*10^-4;
% Resulting max stresses
sigma_r_max=[-5.3537 -3.1726 -2.1168 -1.5973 -1.3396 -1.2113 -1.1473 -1.1153]*10^6;
sigma_theta_max=[1.5044 1.5271 1.5389 1.5448 1.5478 1.5493 1.5501 1.5505]*10^8;

% Plot convergence
figure(4)
plot(nel_conv,u_max,'go-','LineWidth',2)
title('Convergence of radial displacement')
xlabel('Number of elements')
ylabel('max({u_h})')
ylim([3.1985 3.1995]*10^-4)

figure(5)
plot(nel_conv,sigma_r_max,'go-','LineWidth',2)
title('Convergence of radial stress')
xlabel('Number of elements')
ylabel('max({\sigma_r})')

figure(6)
plot(nel_conv,sigma_theta_max,'go-','LineWidth',2)
title('Convergence of tangential stress')
xlabel('Number of elements')
ylabel('max({\sigma_{theta}})')





