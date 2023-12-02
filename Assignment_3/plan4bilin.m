  function [ Ke, fe ] = plan4bilin( ex, ey, ep, D, eq)
% [ Ke, fe] = plan4bilin( ex, ey, ep, D, eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute the stiffness matrix and element external force vector
%  for a bilinear plane stress or plane strain element including 
%  influence of out-of-plane stress (or strain)
%
% INPUT:  ex = [ x1 x2 x3 x4 ]         element nodal x-coordinates
%         ey = [ y1 y2 y3 y4 ]         element nodal y-coordinates
%
%         ep = [ t ngp]                t: thickness
%                                      ngp: number of Gauss points in each
%                                           direction ( ksi and eta)
%                                           (ngp = 1 or 2 or 3)
%
%         D(3,3)                       constitutive matrix for 2D
%                                      elasticity (plane stress or plane
%                                      strain)
%
%         eq = [ bx;                   bx: body force x-dir
%                by ]                  by: body force y-dir
%
%--------------------------------------------------------------
% OUTPUT: Ke : element stiffness matrix (8 x 8)
%         fe : equivalent nodal forces (8 x 1)
%--------------------------------------------------------------
% Developed by Peter Möller Dept. Applied Mechanics, Chalmers
%
% MODIFIED for VSM167 FEM BASICS by Dimosthenis Floros 20151201
%
% MODIFIED for VSM167 FEM BASICS by Martin Fagerström 20211201
%
%--------------------------------------------------------------
%

t     = ep(1);
ngp   = ep(2)^2; % Total gauss points = ( NoGaussPoits per direction )^2

%  Initialize Ke and fe with zeros for all of their elements

Ke = zeros(..........
fe = zeros(..........

%  Tolerance for the Jacobian determinant

minDetJ = 1.e-16;

% Determine constitutive matrix D for plane strain or plane stress

%  Set the appropriate Gauss integration scheme for 1, 2 or 3 Gauss points
%  in each direction (ksi, eta). (or else 1, 4 or 9 Gauss points in total)

if ngp == 1 % 1x 1 integration
 
 intWeight   = .........;
 GaussPoints = .........;

elseif ngp == 4 % 2 x 2 integration
 
 intWeight = ...........;
             
      
 GaussPoints = .........;
                
         
elseif ngp == 9 % 3 x 3 integration
 
 intWeight = ...........;
              
      
 GaussPoints = .........;
               
         
else
 
 error('Only 1,2 or 3 Gauss Points in each direction apply')

end

%  Loop over all integration points to compute Ke and fe 

for gpIndex = 1:ngp
 
 xsi       = ........; 
 weightXsi = ........;
 eta       = ........;
 weightEta = ........;
 
% Compute the element shape functions Ne (use xsi and eta from above)

 Ne = .....;

% Compute derivatives (with respect to xsi and eta) of the
% shape functions at coordinate (xsi,eta). Since the element is
% isoparametic, these are also the derivatives of the basis functions.

 .........

%  Use shape function derivatives and element vertex coordinates 
%  to establish the Jacobian matrix.

 ........

%  Compute the determinant of the Jacobian and check that it is OK

 detJ = det(......);

 if ( detJ < minDetJ )

  fprintf( 1, 'Bad element geometry in function plan4bilin: detJ = %0.5g\n', detJ );
  return;

 end

% Determinant seems OK - invert the transpose of the Jacobian

 .........

%   Compute derivatives with respect to x and y, of all basis functions,

 .....

% Use the derivatives of the shape functions to compute the element
% B-matrix, Be

 Be = ......;

% Compute the contribution to element stiffness matrix and volume load vector 
% from current Gauss point
% (check for plane strain or plane stress again!)     

 Ke = Ke + .....;

 fe = fe + .....;   

end  

end

%--------------------------end--------------------------------

