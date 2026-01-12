% form_forceF
%
function[matrF] = form_forceF(matrF)
%
% This function forms the global force vector
%
global nm mdof
global nf load
%
for i = 1:nm
    for j = 1:mdof
        if nf(i,j) ~= 0
            matrF(nf(i,j)) = load(i,j);
        end
    end
end
%
%%%%%%%%%%%%%%%% end of function form_forceF  %%%%%%%%%%%%%%