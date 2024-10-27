function v=Warp_Image(f,vecX,vecY,interp_method)
%===============================================================================
% function v=Warp_Image(f,vecX,vecY,interp_method)
%
% This function does the morphing of the image via the given vector field
% accordingly to the provided interpolation method.
%
% Inputs:
%   f: image to warp
%   vecX : horizontal component of the deformation field
%   vecY : vertical component of the deformation field
%   interp_method: interpolation method to be used. It uses the methods
%   available in the interp2 function:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - bilinear interpolation
%       'spline'  - spline interpolation
%       'cubic'   - bicubic interpolation as long as the data is
%                      uniformly spaced, otherwise the same as 'spline'
%       'makima'  - modified Akima cubic interpolation
%
% Output:
%   v: Warped image
%
% Author: Jerome Gilles
% Institution: San Diego State University - Dept of Mathematics & Statistics
% Version: 1.0
%===============================================================================

[X,Y]=meshgrid(1:size(f,2),1:size(f,1));

%Set the coordinates where to map
XX=X+vecX(:,:);
YY=Y+vecY(:,:);

%check if the coordinates go outside the image domain
XX(XX<1)=1;
XX(XX>size(f,2))=size(f,2);
YY(YY<1)=1;
YY(YY>size(f,1))=size(f,1);

%maps the image using interpolation
v(:,:)=interp2(f(:,:),XX,YY,interp_method);    
