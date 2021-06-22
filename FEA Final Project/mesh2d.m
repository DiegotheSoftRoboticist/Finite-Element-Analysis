%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [coord, connect] = mesh2d(mesh_type, xmin, ymin, xmax, ymax, 
%                                       ndiv_x, ndiv_y)
%
% Purpose
% =======
% Finite element mesh for rectangular domains with option to
% create either rectangular or triangular (rectangle is sub-
% divided into two) elements. Output plot is also generated to 
% show the mesh with the nodal and element numbering
%
% Function Syntax
% ===============
% function [coord, connect] = 
% mesh2d(mesh_type, xmin, ymin, xmax, ymax, ndiv_x, ndiv_y)
%
% Inputs
% ======                               
% mesh_type   = 'TRI3'  for 3-node triangles and 
%               `QUAD4' for quadrilaterals (rectangles)
%                                               _
% xmin        = x-coordinate of the lower left   | coordinates of
% ymin        = y-coordinate of the lower left   | the domain
% xmax        = x-coordinate of the upper right  | (rectangle)
% ymax        = y-coordinate of the upper right _|
% ndiv_x      = number of sub-divisions in the x-direction
% ndiv_y      = number of sub-divisions in the y-direction
%
% Outputs
% =======        
% coord   = 2-d array which stores the coordinates of the nodes
%           coord(I,1) and coord(I,2) : x- and y-coordinate
%           of node I                      
% connect = element connectivity 
%           this is used to assemble the GSM in the main code
%           row i of this matrix contains the nodes in element i 
%           (in the CCW direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  function [coord, connect] = mesh2d(mesh_type, xmin, ymin, xmax, ymax, ndiv_x, ndiv_y)

  if ~(strcmp(mesh_type, 'QUAD4'))
    error('Mesh generator only coded for ''QUAD4'' element');
    return;
  end

  nodes_x = ndiv_x + 1; 
  nodes_y = ndiv_y + 1; 
  numnd = nodes_x*nodes_y; % numnd is the number of nodes in the mesh

% nodes are equi-spaced in x-direction; hx = x-spacing
  x = linspace(xmin, xmax, nodes_x); hx = x(2) - x(1); 
% nodes are equi-spaced in y-direction; hy = y-spacing
  y = linspace(ymin, ymax, nodes_y); hy = y(2) - y(1);

  for i = 1:nodes_x
    indices = i:nodes_x:(nodes_x*(nodes_y-1) + i);
    coord(indices,1) = x(i).*ones(nodes_y,1); % same x-coordinate along line x = constant
    coord(indices,2) = y(1:nodes_y)';         % y-coordinate from array y()
  end

  
  
%
% code the connectivity 
%



%-------------------------------------------------------------
% Mesh plotting
%-------------------------------------------------------------

  hold on;
  axis([xmin xmax ymin ymax]); % set the axes for plotting
  xlabel('X','fontsize', 12);  % label for x-axis
  ylabel('Y','fontsize', 12);  % label for y-axis

  for j = 1:nodes_y
    ypoints = y(j)*ones(1,nodes_x);
% element edges in `blue' color 
    plot(x,ypoints,'b-');
  end
  for i = 1:nodes_x
    ypoints = linspace(ymin, ymax, nodes_x);
% element edges in `blue' color
    plot(x(i)*ones(1,nodes_x), ypoints,'b-'); 
  end

  if (strcmp(mesh_type,'TRI3'))
%
% code the additional edges for the TRI3 elements and then plot them
%
  end

%
% set the font size for the node and element numbering
%
  numel = length(connect(:,1));
  if (numel <= 50) fsize = 10;
  else             fsize = 8;
  end

%
% node numbering
%
  offset_x = 0.06*hx; offset_y = 0.08*hy;
  for i = 1:nodes_x*nodes_y
    idstring = sprintf('%d',i);
    text(coord(i,1)+offset_x, coord(i,2)+offset_y, idstring,'color',[1,0,0],'fontsize',fsize);
  end

% element numbering
%
  if (strcmp(mesh_type,'QUAD4'))
%
% element numbering for QUAD4
%
    for e = 1:length(connect(:,1))
      xc = coord(connect(e,1),1) + hx/2; % element number is plotted in the
      yc = coord(connect(e,1),2) + hy/2; % center of the rectangle
      idstring = sprintf('%d',e);
      text(xc, yc, idstring, 'color',[0,0,0],'fontsize',fsize);
    end
  else 
%
% element numbering for TRI3
%

%
%  code the element numbering for TRI3
%
  end

  hold off;
