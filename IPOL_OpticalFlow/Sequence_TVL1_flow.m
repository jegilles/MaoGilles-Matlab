function flow=Sequence_TVL1_flow(seq1,tau,lambda,theta,zfactor,epsilon,nwarps,nscales)

flow=zeros(size(seq1,1),size(seq1,2),2,size(seq1,3)-1);

for t=1:size(seq1,3)-1
   [flow(:,:,1,t),flow(:,:,2,t)]=IPOL_TVL1_OpticalFlow(seq1(:,:,t),seq1(:,:,t+1),tau,lambda,theta,zfactor,epsilon,nwarps,nscales);
end