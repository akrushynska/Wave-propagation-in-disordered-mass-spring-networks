%   spring_g.m
%
function[g] = spring_g(i)
%
% This function forms the steering vector for element i
%
global connec nf
%
% retrieve the nodes of element i
%
node1 = connec(i,1);
node2 = connec(i,2);
%
% form the steering vector from element's DOFs
%
g = [nf(node1,1); nf(node1,2); nf(node2,1); nf(node2,2)];
%
%%%%%%%%%%%%%%%% end function spring_g %%%%%%%%%%%%%%
