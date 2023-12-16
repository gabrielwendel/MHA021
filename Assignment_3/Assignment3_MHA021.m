%% Assignment 3 MHA021
% Group 33
% Nils Helgesson & Gabriel Wendel
close all;
clear all;
clc;

%% Indata
E=210e9; %[Pa]
nu=0.3;
rho=8000; %[kg/m3]
g=9.82; %[m/s^2]
b=[0; -rho*g]; %load, self weight
p_max=200e6; %[Pa]
% p_max=0; %For control calc

%Dimensions of half beam
W=0.5;
H=0.3;
% W=10*H; %For control calc
t=0.05;

% Meshdata
mesh_a=40e-3;
nelH=20;
nelW=20;
nela=20;

%Gauss points
ngp=1;



%% Generate mesh
% INPUT: H      : domian height
%        W      : domain width
%        a      : fraction of the width along the top edge
%        nelH   : number of element edges in vertical direction
%        nelW   : number of element edges in horizontal direction
%        nela   : number of element edges along the distance a (nela<nelvW)
%
% OUTPUT:  xy   : matrix with 2 columns; ith row contains the (x,y)-
%                 -coordinate of the ith node (CALFEM name Coord)
%          Dof  : matrix (2 columns) vector with degrees of freedom in
%                 each node; (row i simply has the the two values 2*i-1
%                 and 2*i), since there are two dofs in each node, i.e.
%                 the degrees of freedom in node i, are 2*i-1 and 2*i
%          Edof : matrix with 9 columns; 1st column is the element 
%                 number; columns 2-9 in row i, contains the degrees
%                 of freedom, in counter clock-wise order, for element i.
%  SEE THE CALFEM MANUAL, PAGES 6.2-2 TO 6.2-3, FOR A MORE DETAILED
%  DESCRIPTION OF THE ABOVE OUTPUT VARIABLES.
%          nodL  : list of nodes along left edge
%          nodR  : list of nodes along right edge
%          noda  : list of nodes along part "a" of upper edge



[xy,Dof,Edof,nodL,nodR,noda] = pmesh(H,W,mesh_a,nelH,nelW,nela);

nel=size(Edof,1);
nnodes=size(Dof,1);

%Process mesh data
ex=xy(:,1);
ey=xy(:,2);
Ex=zeros(nel,4);
Ey=zeros(nel,4);
for i=1:nel
    E1=find(ismember(Dof, Edof(i,2:3),'rows'));
    E2=find(ismember(Dof, Edof(i,4:5),'rows'));
    E3=find(ismember(Dof, Edof(i,6:7),'rows'));
    E4=find(ismember(Dof, Edof(i,8:9),'rows'));
    Ex(i,:)=[ex(E1) ex(E2) ex(E3) ex(E4)];
    Ey(i,:)=[ey(E1) ey(E2) ey(E3) ey(E4)];   
end

eldraw2(Ex,Ey)
%% Constitutive matrix D
ptype=1; % plane stress
% Determine constitutive matrix using Hooke's generalized law
D=hooke(1,E,nu);

%% Assemble fb

%along r5
GP=0; %middle of interval
iW=2;
Ne_xsi=[0.5 0 0.5 0; 0 0.5 0 0.5]; %middle of 1D element
fb=zeros(nnodes*2,1);
for i=1:length(noda)-1
    x1=ex(noda(i));
    x2=ex(noda(i+1));
    L_e=x2-x1;
    x_middle=(x1+x2)/2;
    ty=-p_max*sqrt(1-(x_middle/mesh_a)^2);
    fbe=Ne_xsi'*[0 ty]'*t*L_e;
    fb(Dof(noda(i),1))=fbe(1); %x dof for node 1
    fb(Dof(noda(i),2))=fbe(2); %y dof for node 1
    fb(Dof(noda(i+1),1))=fbe(1); %x dof for node 2
    fb(Dof(noda(i+1),2))=fbe(2); %y dof for node 2   
end

%% Assemble K and fl
K=zeros(nnodes*2,nnodes*2);
fl=zeros(nnodes*2,1);
ep=[t ngp];
for i=1:nel
    [Ke fe]=plan4bilin(Ex(i,:),Ey(i,:), ep, D, b);

    [K, fl]=assem(Edof(i,:),K,Ke,fl,fe);
end

% Boundary conditions

bc=[Dof(nodR,1) zeros(length(nodR),1);
    Dof(nodR,2) zeros(length(nodR),1);
    Dof(nodL,1) zeros(length(nodR),1)];
[a , r]=solveq(K , fl+fb , bc ) ;


% Plots
ed = extract ( Edof , a ) ; 
figure (1)
eldraw2( Ex , Ey , [1 4 0])
hold on
eldisp2( Ex , Ey , ed , [1 2 1] , 300) ; % scale factor =300
% eldisp2( Ex , Ey , ed , [1 2 1], 1) ; %scale factor =0
xlabel('x[m]')
ylabel('y[m]')
legend('Undeformed')
title('Deformation')

%% Calculate stresses

nel = size(Edof,1); % number of elements

ep = [ptype,t,ngp];

% Stress vector for element center
sigma = zeros(nel, 3);

for i=1:nel
    % # columns in es follows size of D (3x3) and # rows = n^2, n =
    % integration points
    [es,et,eci] = plani4s(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
    % Store stresses for each element
    sigma(i,:) = es';
end

%% Effective stress (von Mises)

sigma_vm = zeros(1,nel);

for i=1:nel
    % calculate effective stress according to von Mises
    sigma_vm(i) = sqrt(sigma(i,1)^2 + sigma(i,2)^2 - sigma(i,1)*sigma(i,2) + 3*sigma(i,3)^2);
end

figure(2)
xlabel('x-coordinate [m]');
ylabel('y-coordinate [m]');
title('Effective stress (von Mises) [Pa] for each element');

hold on

for i=1:nel
    fill(Ex(i,:),Ey(i,:),sigma_vm(i))
end
colorbar;
colormap('jet')

maxSigma=max(sigma_vm);

