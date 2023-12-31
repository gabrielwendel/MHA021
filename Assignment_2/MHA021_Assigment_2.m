%% MHA021 FEM - Assigment 2
% Group 33
% Gabriel Wendel, Nils Helgesson
close all;
clear all;
clc;

%% Task 1

load('Mesh_data.mat')

%indata
h1=0.35;
h2=0.15;
h3=0.018;
h4=0.2;

T_in=25;
T_out=-10;
T=[T_out T_in]

alpha=10;
L1=1.2;
L2=0.6;
H=1.5;

k1=1.5; %concrete
k2=0.035; %insulation
k3=0.17; %plaster
k=[k1 k2 k3]

thickness=1;

K=zeros(NoDofs);

% Assemble K elementwise
for element = 1:NoElem
    D=k(matrlIndex(element)).*eye(2);
    
    Ke = flw2te(Ex(element,:),Ey(element,:),thickness,D);
    K = assem(Edof(element,:),K,Ke); 
end

%Assemble K with Kc and fb

NoBoundary=length(boundaryEdof)
Kc=zeros(NoDofs);
fb=zeros(NoDofs,1);

for element= 1:NoBoundary
    ex=boundaryEx(element,:);
    ey=boundaryEy(element,:);
    Tamb=T(boundaryMaterial( element, 2 ));
    [Kce, fce] = convecte(ex, ey, alpha, thickness, Tamb);
    [Kc, fb]=assem(boundaryEdof(element,:),Kc,Kce,fb,fce);
end

a=solveq(K+Kc,fb)

% plot
figure(1)
ed=extract(Edof,a);
fill(Ex',Ey',ed')
colormap parula
colorbar

axis equal
axis off

%Linear transmittance
Q_out=0
outside_boundary_index=find(boundaryMaterial(:,2)==1)
for i=1:length(outside_boundary_index)
    ex=boundaryEx(i,:);
    ey=boundaryEy(i,:);
    i_index=boundaryEdof(i,2);
    j_index=boundaryEdof(i,3);
    T_i=a(i_index);
    T_j=a(j_index);
    L_ex=ex(2)-ex(1);
    L_ey=ey(2)-ey(1);
    L_e=sqrt(L_ex^2+L_ey^2);
    Tamb=T(boundaryMaterial(i, 2 ));
    Q_e=alpha*L_e*thickness*((T_i+T_j)/2-Tamb);
    Q_out=Q_out+Q_e;
end
Q_out=2*Q_out;
%analytic solution
R_conv=1/(alpha);
R_1=h1/(k1);
R_2=h2/(k2);
R_3=h3/(k3);
R_ekv=2*R_conv+R_1+R_2+R_3;
q_analytic=(T_in-T_out)/R_ekv

psi=(abs(Q_out)-q_analytic*H)/(T_in-T_out)

%% Task 3

% Lengths [m]
L1 = 4;
L2 = 2.5;
L3 = 2;

% Angle
alpha = 30;     % [deg]

% Weight of gondola plus passengers
m = 3000;       % [kg]     
g = 9.82;       % [m/s^2]

% Point force [N]
P1 = 4000;
P2 = 6000;
P_m = m*g;

% Load [Nm^-1]
q_0 = 6000;
q_1 = 3000;

% Young's modulus
E = 210e9;      % [Pa]

% Beam areas [HEA100 HEA120] [m^2]
A = [21.24*10^-4 25.34*10^-4];

% Beam moment of inertia [HEA100 HEA120] [kg*m^2]
I = [349.2*10^-8 606.2*10^-8];

% Number of beams
num_HEA100_beams = 4;
num_HEA120_beams = 3;

% Matrices containing beam properties
HEA100_prop = repmat([E A(1) I(1)], num_HEA100_beams, 1);
HEA120_prop = repmat([E A(2) I(2)], num_HEA120_beams, 1);

% ep = [E A I]
Ep = zeros(7,3);

% Beam 1,2,3 and 4 are HEA100 beams
Ep(1:4, :) = HEA100_prop;
% Beam 5,6 and 7 are HEA120 beams
Ep(5:end,:) = HEA120_prop;

% Load vector q_e rows = number of beams
q_e = zeros(7,2);
% Beam 5,6 and 7 experience loads
q_e(5,:) = [0,-q_0];
q_e(6:7,:) = [0,-q_1; 0,-q_1];


% Elements dofs
Edof = [1 1:6
        2 4:9
        3 10:15
        4 13:18
        5 4 5 6 13 14 15
        6 7 8 9 16 17 18
        7 16:21];

Dof = [1:3
       4:6
       7:9
       10:12
       13:15
       16:18
       19:21];

% Coordinates
Coords = [0 0; -L1*sind(alpha) L1*cosd(alpha);0 2*(L1*cosd(alpha));L2 0;...
    L2+L1*sind(alpha) L1*cosd(alpha);L2 2*(L1*cosd(alpha));L2+L3 2*(L1*cosd(alpha))];

% Element coordinates Ex, Ey
[Ex, Ey] = coordxtr(Edof, Coords, Dof, 2);

% Illustrate the construction
figure(1)
eldraw2(Ex, Ey, [1 1 1], Edof(:,1))
hold on

% Initialize load vector and stiffness matrix
f = zeros(21,1);
K = zeros(21);

% Compute element stiffness matrix and element load vector for a two 
% dimensional beam element by using CALFEM function "beam2e"
for i=1:length(Edof)
    % Loop over all elements and construct stiffness matrix and load vector
    [Ke,fe] = beam2e(Ex(i,:), Ey(i,:), Ep(i,:), q_e(i,:));
    % Assembly
    [K,f] = assem(Edof(i,:), K, Ke, f, fe);
end

% Add point forces to load vector
f(4) = P2;
f(7) = P1;
f(20) = -P_m;

% Boundary conditions
bc = [1 0; 2 0; 3 0; 10 0; 11 0; 12 0];

% Solve for displacement, a
a = solveq(K, f, bc);

% Maximum displacement
max_disp = max(abs(a));

% Extract all degrees of freedom for each element
Ed = extract_dofs(Edof, a);

sfac = scalfact2(Ex, Ey, Ed, 0.1);

plotpar = [2 1 0];

figure(1)
hold on
eldisp2(Ex, Ey, Ed, plotpar, sfac)
xlabel('x [m]')
ylabel('y [m]')
title('Gondola displacement (ratio = 0.1)')



n = 20; % number of evaluation points along the beam

% Compute sectional forces along the element [N V M]
% column 1 is normal forces, column 2 is th shear force and column 3 is 
% the bending moment 
for i=1:length(Edof)
    SectionalForces(i).es = beam2s(Ex(i,:), Ey(i,:), Ep(i,:), Ed(i,:), q_e(i,:), n);
end

% Scaling factors
sfacs = [5e-5, 3e-5, 2e-5];  

% Titles for the three plots
titles = {'Normal force', 'Shear force', 'Bending moment'};

% plot sectional forces
for j=1:3
    figure(j+1)
    for i=1:length(Edof)
        eldia2(Ex(i, :), Ey(i, :), SectionalForces(i).es(:,j), [2 1], sfacs(j));
        hold on
    end
    xlabel('x [m]')
    ylabel('y [m]')
    title(titles{j})
end

% Maximum bending moment
max_M_beam = zeros(length(Edof),1);
% Loop over structure, max_M_beam = max bending moment for every beam
for i=1:length(Edof)
    max_M_beam(i) = max(abs(SectionalForces(i).es(:,3)));
end

% Max bending moment in structure
max_M = max(abs(max_M_beam));

% z = h/2 for HEA100 and HEA120 beam [m^3]
z = [4.8e-2 5.7e-2];

% Initialize stress matrices, rows = number of beams, columns = number of
% evalution points
sig_upper = zeros(length(Edof), n);
sig_lower = zeros(length(Edof), n);

% Start with HEA100 beam
k = 1;

% Calculate stress on every beam
for i = 1:length(SectionalForces)
    % At every evaluation point on the beam
    % If beam number > 4 => HEA120 beam
    if i > 4
        k = 2;
    end
    for j = 1:n
        sig_upper(i, j) = SectionalForces(i).es(j, 3) / I(k) * z(k);
        sig_lower(i, j) = -SectionalForces(i).es(j, 3) / I(k) * z(k);
    end
end

% Determine max stress and where it occurs
[sig_upper_max, pos_upper] = max(abs(sig_upper(:)));
[sig_lower_max, pos_lower] = max(abs(sig_lower(:)));

