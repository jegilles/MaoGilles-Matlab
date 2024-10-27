function [flv,flh]=IPOL_HornSchunck_OpticalFlow(im1,im2,alpha,zfactor,tol,nwarps,nscales,maxiter)

%******************************************************************************************************
% function [flv,flh]=IPOL_HornSchunck_OpticalFlow(im1,im2,alpha,zfactor,tol,nwarps,nscales,maxiter) 
%
% This function computes the optical flow between two images via the pyramidal Horn-Schunck algorithm.
% It corresponds to the flow which maps im2 to im1.
% It is a Matlab adaptation via mex files for the C implementation from the
% Image Processing OnLine (IPOL) website. Only the needed files were kept
% from the original IPOL archive (see below for a link to the original
% algorithm) and only the float format was replaced by a double format in
% the C code.
%
% Inputs:
%  -im1: first image (must be in double format)
%  -im2: second image (must be in double format)
%  -alpha: regularization parameter (20.0)
%  -zfactor: zoom factor (0.5)
%  -tol: tolerance (0.001)
%  -nwarps: number of warpings (10)
%  -nscales: number of scales in the pyramid (10)
%  -maxiter: maximum number of iterations (150)
%
% Outputs:
%  -flv: vertical component of the optical flow
%  -flh: horizontal component of the optical flow
%
%  Source: IPOL (https://doi.org/10.5201/ipol.2013.20) + J.Gilles
%  Date: 06/20/2017
%******************************************************************************************************