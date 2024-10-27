function [vecX,vecY]=Staboptflow(reference,frames,optflow_method)
%===============================================================================
% function [vecX,vecY]=Staboptflow(reference,frames,optflow_method)
%
% Extract the deformation fields via optical flow which map
% the reference image to each frame
%
% Inputs:
% - reference: reference image
% - frames: 3D matrix sequence of frames (3rd dimension is time)
% - optflow_method: choice of which optical flow method to use, the choices are
%                       - 'lk' : Lucas-Kanade (you ned Piotr toolbox to run this
%                       option)
%                       - 'tvl1' : TV-L1 model
%                       - 'hs' : Horn-Schunck model
%
% Outputs:
%   vecX : sequence of horizontal displacements
%   vecY : sequence of vertical displacements
%
% Author: Yu Mao and Jerome Gilles
% Version 2.0
%===============================================================================

S1=size(frames,1);
S2=size(frames,2);
S3=size(frames,3);

vecX=zeros(S1,S2,S3);
vecY=zeros(S1,S2,S3);

eps=1e-8;

for k=1:S3
    
    switch lower(optflow_method)
        case 'lk'
        %Lucas-Kanade optical-flow
        [uvX,uvY,~] = optFlowLk(frames(:,:,k),reference,[],4,1.2,3e-6,0);
        
        case 'tvl1'
        %TV-L1 optical-flow
        [uvY,uvX]=IPOL_TVL1_OpticalFlow(frames(:,:,k),reference, ...
            0.25,0.15,0.3,0.5,0.01,5,10);
        
        case 'hs'
        %Horn-Schunck optical-flow
        [uvY,uvX]=IPOL_HornSchunck_OpticalFlow(frames(:,:,k),reference, ...
            20.0,0.5,0.001,10,5,10);
    end
       
    % Take care of the boundary such that no vector points out of the region.
    [X,Y]=meshgrid(1:S2,1:S1);
    
    OUTX=(X+uvX>S2-eps);
    uvX(OUTX)=S2-X(OUTX)-eps;
    OUTX=(X+uvX<1+eps);
    uvX(OUTX)=1-X(OUTX)+eps;
    vecX(:,:,k)=uvX;
    
    OUTY=(Y+uvY>S1-eps);
    uvY(OUTY)=S1-Y(OUTY)-eps;
    OUTY=(Y+uvY<1+eps);
    uvY(OUTY)=1-Y(OUTY)+eps;
    vecY(:,:,k)=uvY;
end