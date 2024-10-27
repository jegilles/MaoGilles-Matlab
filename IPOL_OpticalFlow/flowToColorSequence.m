function cseq=flowToColorSequence(u,v)

%===============================================================================
% function cseq=flowToColorSequence(u,v)
%
% This function return a colorized version of the optical flow of a sequence (it
% uses the flowToColor function).
%
% Inputs:
%   u: 3D matrix containing the horizontal flow (the third dimension corresponds
%   to time
%   v: 3D matrix containing the vertical flow (the third dimension corresponds
%   to time
%
% Output:
%   cseq: 4D matrix containing the colorized optical flow (the fourth dimension
%   corresponds to time
%
% Author: Jerome Gilles
% Institution: San Diego State University - Dept of Mathematics & Statistics
% Version: 1.0
%===============================================================================

cseq=zeros(size(u,1),size(u,2),3,size(u,3));
for k=1:size(u,3)
    cseq(:,:,:,k)=flowToColor(u(:,:,k),v(:,:,k));
end