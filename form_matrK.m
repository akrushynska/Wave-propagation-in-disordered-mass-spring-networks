% form_matrK
%
function[matrK] = form_matrK(matrK, kg,g)
%
% This function assembles the global stiffness materix
%
global spdof
%
for i = 1:spdof
    if g(i) ~= 0
        for j = 1:spdof
            if g(j) ~= 0 
                matrK(g(i),g(j)) = matrK(g(i),g(j)) + kg(i,j);
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%% end of functions form_matrK %%%%%%%%%%%%%%%%