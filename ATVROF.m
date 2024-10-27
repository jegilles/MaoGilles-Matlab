function u=ATV_ROF(f,mu,lambda,Niter)

%===================================================================
% function u=ATV_ROF(f,mu,lambda,Niter)
%
% Anisotropic ROF denoising
% Version:
% -v1.2 - 06/20/2013
% -v1.0 - 04/04/2011
%
% This function performs the minimization of the Anisotropic TV - L2 
% Rudin-Osher-Fatemi model for images by Split Bregman Iterations
%
% f = noisy image
% mu = regularization parameter
% lambda = regularization for the Bregman variable
% Niter = maximum number of iterations
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%===================================================================

[M,N]=size(f);

f=double(f);
dx=zeros(M,N);
dy=zeros(M,N);
bx=zeros(M,N);
by=zeros(M,N);
u=f;
Z=zeros(M,N);

Mask=zeros(M,N);
Mask(1,1) = 1;
FMask=fft2(Mask);

%Fourier Laplacian mask initialization
D = zeros(M,N);
D([end,1,2],[end,1,2]) = [0,1,0;1,-4,1;0,1,0];
FD=fft2(D);

%Fourier constant initialization
FW=((mu/lambda)*abs(FMask).^2-real(FD)).^-1;
FF=(mu/lambda)*conj(FMask).*fft2(f);

K=1;
err=norm(f(:),2);
tol=1e-3*norm(f(:),2);
while ((err>tol) && (K<Niter)),
    K=K+1;
    tx=dx-bx;
    ty=dy-by;
    up=u;
    %Update u
    u=real(ifft2(FW.*(FF-fft2(tx-tx(:,[1,1:N-1])+ty-ty([1,1:M-1],:)))));
    ux=u-u(:,[1,1:N-1]);
    uy=u-u([1,1:M-1],:);
        
    %Update dx
    tmpx=ux+bx;
    dx=sign(tmpx).*max(Z,abs(tmpx)-1/lambda);
    
    %Update dy
    tmpy=uy+by;
    dy=sign(tmpy).*max(Z,abs(tmpy)-1/lambda);
   
    %Update bx and by 
    bx=tmpx-dx;
    by=tmpy-dy;
    
    err=sum(sum((up-u).^2));
end