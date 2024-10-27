# MaoGilles-Matlab
 Mao-Gilles turbulence stabilization algorithm

 This code implement the stabilization algorithm described in Y.Mao, J.Gilles, "Non rigid geometric distortions 
 correction - Application to Atmospheric Turbulence Stabilization", Inverse Problems and Imaging Journal, Vol.6, 
 No.3, 531-546, Aug. 2012.

 The main functions are:
 - MaoGillesStabilization.m which return a stabilized image from a sequence of distorted images
 - MaoGillesStabilization_Seq.m which return a sequence of stabilized images from a sequence of distorted images

 These functions use some optical flow estimators which are provided by:
 - The IPOL_OpticalFlow folder: these are Matlab wrappers for the C code available at https://www.ipol.im/
 - the Piotr's Toolbox (v2.52 provided as a zip file in this package for compatibility purposes, check their licence)
