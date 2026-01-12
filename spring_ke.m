% spring_ke.m
%
function[ke] = spring_ke(i)
%
% This function forms the element stiffness matrix in local coordinates
%
global geom connec prop
%
% retrieve the nodes of element i
%
node1 = connec(i,1);
node2 = connec(i,2);
%
% Retrieve the node coordinates
%
x1 = geom(node1,1); y1 = geom(node1,2);
x2 = geom(node2,1); y2 = geom(node2,2);
%
% Evaluate length of element i
%
L = sqrt((x2-x1)^2 + (y2-y1)^2);
%
% Estimate the stiffness of element i
%
E = prop*L;
%
% Calculate element stiffness matrix in local coordinates
%
ke = [E  0 -E  0; ...
      0  0  0  0; ...
     -E  0  E  0; ...
      0  0  0  0];
  %
 %%%%%%%%%%%%%%%%%%%% End of function spring_ke.m %%%%%%%%%%%%%%%%%%%%