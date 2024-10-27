function wheel = GenerateColorWheel(s)

%===============================================================================
% function wheel = GenerateColorWheel(s)
%
% This function create the colorwheel where the color and its intensity 
% correspond to direction and magnitude of optical flows.
%
% Input:
%   s: the size of the wheel will be 2s+1
%
% Ouput:
%   wheel: colorwheel (color matrix)
%
% Author: Jerome Gilles
% Institution: San Diego State University - Dept of Mathematics & Statistics
% Version: 1.0
%===============================================================================

u = repmat([-s:s],2*s+1,1);
v = repmat([-s:s]',1,2*s+1);

wheel = flowToColor(u,v);

for i = 1:size(wheel,1)
    for j=1:size(wheel,2)
        if sqrt((i-s-1)^2+(j-s-1)^2) > s+1
            wheel(i,j,:) = 255;
        end
    end
end