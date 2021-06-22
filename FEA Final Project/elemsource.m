%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [fe] = elemsource(e, nel, ndof, nsp, coord, 
%                            connect, gauss, weight)
% Purpose
% =======
% Element force (load) vector for a 4-node quadrilaterial 
% element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fe] = elemsource(e, nel, ndof, nsp, coord, connect, gauss, weight)

%
% initialize element force (load) vector to zero
%


%
% loop over Gauss points in xi and eta-directions
%
  

% shape function derivatives in the parent coordinate system (xi-eta)
    
% derivative of x and y w.r.t. xi and eta (using isoparametric transformation)
    
% shape function
    
%
%     location of Gauss point in the xy-coordinate system
%
%
% assemble element source vector fe
%
