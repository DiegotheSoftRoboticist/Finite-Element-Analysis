%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function post_processing(xmin, xmax, nodes_x, ymin, ymax, 
%                          nodes_y, solution_coeff)
% Purpose
% =======
% Surface plots of FE solution u^h, exact solution u, and
% the error e = u - u^h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post_processing(xmin, xmax, nodes_x, ymin, ymax, nodes_y, solution_coeff)

 x = linspace(xmin, xmax, nodes_x);
 y = linspace(ymin, ymax, nodes_y);
 [X, Y] = meshgrid(x,y);
%
% exact solution 
%
u = 0.5*(1 - X.^2).*(1 - Y.^2); % exact solution

%
% finite element solution at the nodes
%
 uh = zeros(nodes_y, nodes_x); % note the order
 k = 0;
 for i = 1:nodes_y
   for j = 1:nodes_x
     k = k + 1;
     uh(i,j) = solution_coeff(k);
   end
 end
%
% plotting attributes and the surface plot
%
 figure;%(gcf+1); % new figure window (gcf returns the handle to the previous figure)
 rotate3d;      % can rotate the orientation of the figure on the Matlab window
 surf(X,Y,u); legend('u^{exact}'); % plot the exact solution
 disp('Hit any key to continue . . . ');
 pause;

 figure;%(gcf+1); % new window
 rotate3d;      % can rotate the orientation of the figure on the Matlab window
 surf(X,Y,uh); legend('u^h'); % plot the finite element solution
 disp('Hit any key to continue . . . ');
 pause;

%
% point-wise error e = ( u - u^h )
%
 e = u - uh;
 figure;%(gcf+1); % new figure
 rotate3d; 
 surf(X,Y,e); legend('e = ( u - u^h )'); % plot the error
 disp('Hit any key to continue . . . ');
 pause;

%
% maximum absolute point-wise error at the nodes
%
 maxerror = max(max(abs(e),[],1));
 disp('Max absolute error at a node is: '); disp(maxerror);
