%% MHA021 FEM - Assigment 2
close all;
clear all;
clc;

%% Task 1

load('Mesh_dataCoarse.mat')

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










