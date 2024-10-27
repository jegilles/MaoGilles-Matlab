function u=MaoGillesStabilization(f,lambda,Niter,regularizer,optflow_method)

%===============================================================================
% Mao-Gilles Stabilization algorithm for Atmospheric Turbulence
%
% Inputs:
% - f: image sequence stored in 3D array f(x,y,t).
% - lambda: regularization parameter (a default value of 10000 seems to
%           work in most cases)
% - Niter: number of Bregman iterations (a default value of 4 works well)
% - regularizer: choose either 'NLTV' or 'TV'
% - optflow_method: choice of which optical flow method to use, the choices are
%                       - 'lk' : Lucas-Kanade (you need Piotr toolbox to run
%                       this option)
%                       - 'tvl1' : TV-L1 model
%                       - 'hs' : Horn-Schunck model
%
% Author: Yu Mao, Jerome Gilles
% Version: 2.0
%===============================================================================

%We initializa some variables
numframe=size(f,3); % number of frames
u = mean(f,3); % temporal mean as the initial guess
dt=0.5;        %default=0.5
fk=f; % the Bregman variable

% Bregman Iteration
for j=1:Niter
    % Estimate the deformation flows
    [vecX,vecY]=Staboptflow(u,f,optflow_method);
    
    % Operator Splitting
    for iter=1:5
        du=StabPhiAdj((StabPhi(u,vecX,vecY)-fk),vecX,vecY)/numframe;
        u=u-dt*du;

        switch regularizer
            case 'NLTV' %Nonlocal TV optimized C code
                u = NLTV(u,lambda);lambda=lambda+0.1;             
            case 'TV' %Anisotropic TV
                u=ATVROF(u,lambda,100,10); 
        end
    end
    % Bregman Update
    fk=fk+f-StabPhi(u,vecX,vecY);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=StabPhiAdj(v,vecX,vecY)
%===============================================================================
% Adjoint of the morphing operator
% Inputs:
% - v: the sequence to morph
% - vecX: X component of the deformation field
% - vecY: Y component of the deformation field
%===============================================================================

M=zeros(size(vecX));
for k=1:size(vecX,3)
    M(:,:,k)=PushForward(v(:,:,k),vecX(:,:,k),vecY(:,:,k));
end
u=sum(M,3);



function v=StabPhi(u,vecX,vecY)
%===============================================================================
% Morphing of the sequence u via given sequence of vector field (vecX,vecY).
% A linear interpolation is used to find the pixels outside the grid
%
% Inputs:
% - u: the sequence to morph
% - vecX: X component of the deformation field
% - vecY: Y component of the deformation field
%===============================================================================

[X,Y]=meshgrid(1:size(u,2),1:size(u,1));

v=zeros(size(vecX));
for k=1:size(vecX,3)
    v(:,:,k)=interp2(u,X+vecX(:,:,k),Y+vecY(:,:,k),'linear');    
end
