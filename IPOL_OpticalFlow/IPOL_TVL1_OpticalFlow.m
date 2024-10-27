function [flv,flh]=IPOL_TVL1_OpticalFlow(im1,im2,tau,lambda,theta,zfactor,epsilon,nwarps,nscales)

%******************************************************************************************************
% function [flh,flv]=IPOL_TVL1_OpticalFlow(im1,im2,tau,lambda,theta,zfactor,epsilon,nwarps,nscales) 
%
% This function computes the optical flow between two images via the multiscale TV-L1 algorithm.
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
%  -tau: time step (0.25)
%  -lambda:  regularization parameter (0.15)
%  -theta: tighness parameter (0.3)
%  -zfactor: zoom factor (0.5)
%  -epsilon: stopping criterion (0.01)
%  -nwarps: number of warpings (5)
%  -nscales: number of scales in the pyramid (100)
%
% Outputs:
%  -flv: vertical component of the optical flow
%  -flh: horizontal component of the optical flow
%
%  Source: IPOL (https://doi.org/10.5201/ipol.2013.26) + J.Gilles
%  Date: 06/20/2017
%******************************************************************************************************