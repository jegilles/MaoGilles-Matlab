function flow=Sequence_HornSchunck_flow(seq1,alpha,zfactor,tol,nwarps,nscales,maxiter)

flow=zeros(size(seq1,1),size(seq1,2),2,size(seq1,3)-1);

for t=1:size(seq1,3)-1
   [flow(:,:,1,t),flow(:,:,2,t)]=IPOL_HornSchunck_OpticalFlow(seq1(:,:,t),seq1(:,:,t+1),alpha,zfactor,tol,nwarps,nscales,maxiter);
end