function [ Kce, fce ] = convecte( ex, ey, alpha, thick, Tamb )
% function [ Kce, fce ] = convecte( ex, ey, alpha, thick, Tamb )
%-------------------------------------------------------------------
% Purpose: Compute the element stiffness matrix and force vector
%          due to convection for a linear boundary segment
%-------------------------------------------------------------------
% Input: ex      Boundary element nodal x-coordinate, size(ex) = [1, 2]
%        ey      Boundary element nodal y-coordinate, size(ey) = [1, 2]
%        alpha   Coefficient of thermal expansion, [W/(m^2*degC)]
%        thick   Out-of-plane thickness of the continuum
%        Tamb    Ambient temperature at the boundary segment, [degC]
%-------------------------------------------------------------------
% Output: Kce    Boundary element stiffness matrix due to convection
%         fce    Boundary element force vector due to convection
%-------------------------------------------------------------------
% Created by: Martin Fagerström 2020-11-10
%-------------------------------------------------------------------

L_ex=ex(2)-ex(1);
L_ey=ey(2)-ey(1);
L_e=sqrt(L_ex^2+L_ey^2);

syms s

N_e=[1-s/L_e, s/L_e];

Kce=int(N_e'*N_e*alpha,s,0,L_e);
fce=int(N_e'*alpha*Tamb,s,0,L_e);


%Please complete the function file ot kalculate Ke and fce

