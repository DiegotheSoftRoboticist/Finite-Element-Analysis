%-------------------------------------------------------------
% Course ME 554
%
% Instructor: Heidi Feigenbaum
% E-mail    : heidi.feigenbaum@nau.edu
% Date      : Spring 2018
%-------------------------------------------------------------
%
% MODEL 2D Problem (Strong Form)
% ==============================
% Two-dimensional Matlab finite element program to
% solve the following model problem (steady-state
% heat conduction in a plate):
% 
% Div (k Grad u) + f = 0 in \Omega = (XMIN,YMIN) x (XMAX,YMAX), where
% k = k(x,y) is the conductity (plate is assumed to be thermally
% isotropic), and f = f(x, y) is the distributed heat source in the
% domain (no point sources or convection is assumed).
%  
% Boundary conditions:
%
% u = 0 on \Gamma = \partial \Omega (Essential Boundary Conditions)

% FE Approximation
% ================
% Galerkin approximation with isoparametric four-node bilinear 
% quadrilateral elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear            % clear all variables
close all        % close all figures
clc
%
% Input Data
% ==========
% 1.  Domain     = (xmin,ymin) x (xmax,ymax)
% 2.  mesh_type  = 'QUAD4' (4-node quadrilateral element)
% 3.  nel        = number of nodes per element
% 4.  ndof       = number of degrees of freedom per node
% 5.  ndviv_x    = number of sub-divisions in the x-direction 
% 6.  ndviv_y    = number of sub-divisions in the y-direction 

% 7.  ebcdof     = vector consisting of dofs associated with essential
%                  boundary conditions
% 8.  ebcval     = vector consisting of values associated with
%                  the EBCs in ebcdof
% 9.  k(x,y)     = conductivity defined in diffusive.m 
% 10. f(x,y)     = distributed source defined  in source.m

% 11. nsp x nsp  = Gauss quadrature rule used for the
%                  numerical integration

xmin = -1; ymin = -1; xmax = 1; ymax = 1; % problem domain

mesh_type = 'QUAD4';
nel       = 4;
ndof      = 1;
ndiv_x    = 10;
ndiv_y    = 10;

%
% compute nodes in x- and y-directions, total number of nodes and number of elements
%
nodes_x   = ndiv_x + 1; nodes_y = ndiv_y + 1; 
numnd     = nodes_x*nodes_y; % number of nodes
numel     = ndiv_x*ndiv_y;   % number of elements

%
% mesh information: nodal coordinates and connectivity
%
[coord, connect] = mesh2d(mesh_type, xmin, ymin, xmax, ymax, ndiv_x, ndiv_y);
disp('Hit any key to continue . . . ');
pause;

%
% dofs on the boundary of the domain which are associated with the essential 
% boundary condition and their values
%
ebcdof = [1:nodes_x, nodes_x+1:nodes_x:nodes_x*(ndiv_y-1)+1, 2*nodes_x:nodes_x:nodes_x*ndiv_y, nodes_x*ndiv_y+1:nodes_x*nodes_y];
ebcval = zeros(1, length(ebcdof));

%
% initialize external (load) vector and stiffness matrix
%

numeqns = numnd;                  % number of equations
fext = zeros(numeqns,1);          % initialization of fext
bigk = zeros(numeqns, numeqns);   % initialization of K

% Gauss quadrature
% ================
% 2-point rule in xi and eta: Gauss points = [1/sqrt(3), -1/sqrt(3)],
%                             weights      = [1, 1]
%

nsp     = 2;                          % Gauss-rule (nspxnsp-point rule is used)
gauss   = [ -1/sqrt(3) , 1/sqrt(3) ]; % Gauss quadrature points
weight  = [ 1 , 1 ];                  % weights

%
% compute and assemble stiffness matrix K and external force (load) vector f
%

%
% loop over all the elements
%
  for e = 1:numel         
%
% compute element stiffness matrix
%
  ke = elemstiff(e, nel, ndof, nsp, coord, connect, gauss, weight);
%  
% code the assembly into global stiffness matrix (bigk)
%
   


%
% assemble external force (load) vector
%
% compute element force (load) vector
%
    fe = elemsource(e, nel, ndof, nsp, coord, connect, gauss, weight);
    for i = 1:nel
      nd_I = connect(e,i);
      fext(nd_I) = fext(nd_I) + fe(i);
    end
  end

% enforce essential boundary conditions 

  nebc = length(ebcdof);
  for i = 1:nebc
    n = ebcdof(i);
    for j = 1:numeqns
      fext(j) = fext(j) - bigk(j,n)*ebcval(i);
    end
    bigk(n,:) = zeros(1,numeqns);
    bigk(:,n) = zeros(numeqns,1);
    bigk(n,n) = 1.0;
    fext(n) = ebcval(i);
  end

%
% solve the discrete system: Ku = f
%
  solution_coeff = inv(bigk)*fext;

%
% post-processing
%

  post_processing(xmin, xmax, nodes_x, ymin, ymax, nodes_y, solution_coeff)

