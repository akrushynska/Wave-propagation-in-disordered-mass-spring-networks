% matr_C
%
function[C] = matr_C(theta)
%
% This function return a transformation matrix between local and global
% coordinates for a specified angle 'theta'
%
C = [ cos(theta) -sin(theta)      0           0; ...
      sin(theta)  cos(theta)      0           0; ...
          0           0       cos(theta) -sin(theta); ...
          0           0       sin(theta)  cos(theta);];
%
%%%%%%%% end of function matr_C %%%%%%%%%%%%%
