function useq=MaoGillesStabilization_Seq(f,N,lambda,Niter,regularizer,optflow_method)

%===============================================================================
% Apply the Mao-Gilles stabilization on a shifting temporal windows in order
% to get a stabilized sequence
%
% Inputs:
% - f: input sequence stored in 3D array f(x,y,t).
% - N: length of the temporal window
% - lambda: regularization parameter (a default value of 10000 seems to
%           work in most cases)
% - Niter: number of Bregman iterations (a default value of 4 works well)
% - regularizer: choose either 'NLTV' or 'TV'
% - optflow_method: choice of which optical flow method to use, the choices are
%                       - 'lk' : Lucas-Kanade (you ned Piotr toolbox to run this
%                       option)
%                       - 'tvl1' : TV-L1 model
%                       - 'hs' : Horn-Schunck model
%
% Output:
% - useq : stabilized sequence
%
% Author: Yu Mao, Jerome Gilles
% Version: 2.0
%===============================================================================

L=size(f,3);
M=L-N;
useq=zeros(size(f,1),size(f,2),M);
    
for k=0:M-1
    u=MaoGillesStabilization(f(:,:,[1+k:N+k]),lambda,Niter,regularizer, ... 
        optflow_method);
    useq(:,:,k+1) = u;
end