%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ke] = elemstiff(e, nel, ndof, nsp, coord, 
%                           connect, gauss, weight)
% Purpose
% =======
% Element stiffness matrix for a 4-node quadrilaterial 
% element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ke] = elemstiff(e, nel, ndof, nsp, coord, connect, gauss, weight)

%
% initialize element stiffness matrix to zero
%
elementstiff = zeros(nodes_x,nodes_y);
 
%
% loop over Gauss points in xi and eta-directions
%
 for i = 1:e
 
 element(i)=i,
 
 xi(i,1)= -1;
 eta(i,1)=-1;
 
 xi(i,2)= 1;
 eta(i,2)=-1;
   
 xi(i,3)= 1;
 eta(i,3)=1;
 
 xi(i,4)= -1;
 eta(i,4)=1;
 
 end
 
%
% shape function derivatives in the parent coordinate system (xi-eta)

N1 = (1/4)*(1-xi)*(1-eta)
N2 = (1/4)*(1-xi)*(1+eta)
N3 = (1/4)*(1+xi)*(1+eta)
N4 = (1/4)*(1+xi)*(1-eta)


%
% derivative of x and y w.r.t. xi and eta (isoparametric transformation)
%
for i=1:nodes_x*nodes_y
    
dx_dxi(i) = (1/4)*(-x(i)+x(i+1)+x(i+2)-x(i+4))+(nu/4)*(x(i)-x(i+1)+x(i+2)-x(i+4));
dx_deta(i) = (1/4)*(-x(i)-x(i+1)+x(i+2)+x(i+4))+(xi/4)*(x(i)-x(i+1)+x(i+2)-x(i+4));
dy_dxi(i) = (1/4)*(-y(i)+y(i+1)+y(i+2)-y(i+4))+(nu/4)*(y(i)-y(i+1)+y(i+2)-y(i+4));
dy_deta(i) = (1/4)*(-y(i)-y(i+1)+y(i+2)+y(i+4))+(yi/4)*(y(i)-y(i+1)+y(i+2)-y(i+4));

end 

%  
% compute inverse of the Jacobian matrix; the 1/detj term is not included. Since
% a (1/detj) term arises from J^{-T} (inverse transpose of J), one 
% from J^{-1}, and a detj term from the coordinate transformation, only ONE (1/detj)
% remains which is accounted for in the final result of the element stiffness 
% matrix [see Eq. (*) below]
%

for i=1:nodes_x*nodes_y
J = 1/((dx_dxi(i)*dy_deta(i))-(dx_deta(i)*dy_dxi(i)));
end

      % compute derivatives of shape functions in x- and y-coordinate system

%
% shape function
%

%
% location of Gauss point in the xy-coordinate system
%
x = [N1 N2 N3 N4]*[x(i);x(i+1);x(i+3);x(i+4)]
y = [N1 N2 N3 N4]*[y(i);y(i+1);y(i+3);y(i+4)]

%
% assemble stiffness matrix ke 
%

 %for all elements
    
part1 = [dy_deta(i) -dy_dxi(i) 0 0;-dx_deta(i) dx_dxi(i) 0 0;...
    0 0 dy_deta(i) -dy_dxi(i);0 0 -dx_deta(i) dx_dxi(i)]

part2 = [(-1+eta) 0 (1-eta) 0 (1+eta) 0 (-1-eta) 0;(-1+xi) 0 (-1-xi) 0 (1+xi) 0 (1-xi) 0;...
    0 (-1+eta) 0 (1-eta) 0 (1+eta) 0 (-1-eta);0 (-1+xi) 0 (-1-xi) 0 (1+xi) 0 (1-xi)]


B = (1/(4*J))*[1 0 0 0;0 0 0 1;0 1 1 0]*part1*part2

%Two-dimensional Gauss integration formulas can be obtained by combining two one-di¬
%mensional Gauss quadrature formulas as shown below:



