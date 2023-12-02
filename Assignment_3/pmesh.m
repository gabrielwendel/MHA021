  function[xy,Dof,Edof,nodL,nodR,noda] = pmesh(H,W,a,nelH,nelW,nela)
%
% PURPOSE                                      y
%  Mesh generation, by means of transfinite   |
%  interpolation, on an rectangular domain    |
%  of width W and height H; mesh density is   ------_________------> x
%  enriched towards the top left corner.      |  a     W-a   |
%                                             |              |
%  Output variables tailored to suit CALFEM,  |              |H
%  with the assumption that there are two (2) |              |
%  degrees of freedom in each node.           |              |
%  Bilinear elements.                         |______________|
%  Orign of the coordinate system is at the          W
%  top left corner.
%
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
%
%  Created October 6 2014; PW Möller
%

%  Check, and if necessary, adjust discretiztion parameters
if (nelH < 1) nelH = 1; end
if (nelW < 2) nelW = 2; end
if (nela >= nelW) nela = floor((nelW+1)/2); end

%  Number of elements and nodes to generate: nel and nno
nel = nelH*nelW;
nno = (nelH+1)*(nelW+1);

%  Create output variable Dof
Dof = [[1:2:2*nno]' [2:2:2*nno]'];

%  Allocate space for the other output variables
xy = zeros(nno,2);
Edof = zeros(nel,9);

% Insert element numbers in Edof
Edof(:,1) = [1:nel]';

% work space
xy1 = zeros(2,nelH+1); % lower line node coordinates
xy2 = zeros(2,nelW+1); % right line coordinates
xy3 = zeros(2,nelH+1); % upper line node coordinates
xy4 = zeros(2,nelW+1); % left line node coordinates
% 
% corner coordinates for tfi
c1 = [0  0]';
c2 = [0 -H]';
c3 = [W -H]';
c4 = [W  0]';

% Discretize lower and upper lines
% 
xb = 0; xt = W;
xy1(:,1) = c1; xy1(:,end) = c2;
xy3(:,1) = c4; xy3(:,end) = c3;
dxsi = 1/nelH;
xsi = 0;
if nelH>1
    for i = 2:nelH
        xsi = xsi + dxsi;
        xy1(:,i) = [xb -H*xsi^2]';
        xy3(:,i) = [xt -(i-1)*H*dxsi]';
    end
end

%
% discretize right line
yr = -H; xm = (a + W/2)/2;
xy2(:,1) = c2; xy2(:,end) = c3;
xy2(:,nela+1) = [xm yr]';
deta = xm/nela;
if nela>1
    for i = 2:nela
        xy2(:,i) = [(i-1)*deta yr]';
    end
end
deta = (W - xm)/(nelW-nela);
if (nelW-nela)>1
    for i = 2:(nelW-nela)
        xy2(:,nela+i) = [xm+(i-1)*deta yr]';
    end
end

% discretize left line
yr = 0;
xy4(:,1) = c1; xy4(:,nela+1) = [a yr]'; xy4(:,end) = c4;
% line a
dx = 1/nela;
xr = 0;
if (nela>1)
    for i = 2:nela
        xr = xr + dx*a;
        xy4(:,i) = [xr yr]';
    end
end
% x=a to x=W
xr = a;
dx = (W-a)/(nelW-nela);
if (nelW-nela)>1
    for i = 2:(nelW-nela)
        xr = xr + dx;
        xy4(:,nela+i) = [xr yr]';
    end
end
% xy4
% generate nodes       
% ino = node count
ino = 0;
deta = 1/nelW; eta = -deta;
dxsi = 1/nelH;
for ieta = 1:nelW+1
    eta = eta + deta;
    xsi = -dxsi;
    for ixsi = 1:nelH+1
        xsi = xsi + dxsi;
        ino = ino+1;
        xy(ino,:) = (1-eta)*xy1(:,ixsi)'+eta*xy3(:,ixsi)';
        xy(ino,:) = xy(ino,:) + (1-xsi)*xy4(:,ieta)'+xsi*xy2(:,ieta)';
        te = xsi*eta*c3 + xsi*(1-eta)*c2;
        te = te + eta*(1-xsi)*c4 + (1-xsi)*(1-eta)*c1;
        xy(ino,:) = xy(ino,:) - te';
%         if ixsi==1
%             xy(ino,:)
%         end
    end
end
  
% generate elements - 2 dofs in each node
% iel = element count
iel = 0;
% Loop in eta direction
for ieta = 1:nelW
% First node - at upper left of element; subtract 1 from number
   ni = (ieta - 1)*(nelH + 1);

%  Loop in xsi direction
   for ixsi = 1:nelH
      Edof(iel+ixsi,2:5) = [2*(ni+ixsi)-1 2*(ni+ixsi)  2*(ni+ixsi+1)-1 2*(ni+ixsi+1)];
      Edof(iel+ixsi,6:7) = [2*(ni+ixsi+nelH+2)-1 2*(ni+ixsi+nelH+2)];
      Edof(iel+ixsi,8:9) = [2*(ni+ixsi+nelH+1)-1 2*(ni+ixsi+nelH+1)];
   end
   iel = iel + nelH;
end % end element generation

% selected boundary node lists
nodL = 1:(nelH+1);
noda = 1:(nelH+1):(nela*(nelH+1)+1);
nodR = size(xy,1)-nelH:size(xy,1);

 