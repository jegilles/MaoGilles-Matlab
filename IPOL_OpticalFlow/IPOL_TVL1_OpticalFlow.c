/**************************************************************************
 *
 *  This code interface the DualTVL1OpticalFlow function available at IPOL
 *  The original C source code was modified to handle double instead of 
 *  float for an easier compatibility with Matlab
 *
 *************************************************************************/

#include <matrix.h>
#include <mex.h>  
#include "tvl1flow_lib.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 //declare variables
    mxArray *in_im1, *in_im2, *out_flv, *out_flh,*in_tau,*in_lambda,*in_theta,*in_zfactor,*in_epsilon,*in_nwarps,*in_nscales;
    const mwSize *dims;   

    double *im1,*im2,*flv,*flh,tau,lambda,theta,zfactor,epsilon,N;
    int dimx,dimy,nwarps,nscales;
    
 //associate inputs
    in_im1 = mxDuplicateArray(prhs[0]);
    in_im2 = mxDuplicateArray(prhs[1]);
    in_tau = mxDuplicateArray(prhs[2]);
    in_lambda = mxDuplicateArray(prhs[3]);
    in_theta = mxDuplicateArray(prhs[4]);
    in_zfactor = mxDuplicateArray(prhs[5]);
    in_epsilon = mxDuplicateArray(prhs[6]);
    in_nwarps = mxDuplicateArray(prhs[7]);
    in_nscales = mxDuplicateArray(prhs[8]);

//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    dims = mxGetDimensions(prhs[1]);
    if ((dimx != (int)dims[1])||(dimy != (int)dims[0])) mexErrMsgTxt("The two input images must be of the same size!");

//associate outputs
    out_flv = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    out_flh = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers
    im1 = mxGetPr(in_im1);
    im2 = mxGetPr(in_im2);
    flv = mxGetPr(out_flv); //vertical component of the flow
    flh = mxGetPr(out_flh); //horizontal component of the flow

//associate scalars
    tau = mxGetScalar(in_tau);
    lambda = mxGetScalar(in_lambda);
    theta = mxGetScalar(in_theta);
    zfactor = mxGetScalar(in_zfactor);
    epsilon = mxGetScalar(in_epsilon);
    nwarps = (int)mxGetScalar(in_nwarps);
    nscales = (int)mxGetScalar(in_nscales);
    
 //run the TV-L1 algorithm
    N = 1+log(sqrt(dimx*dimx+dimy*dimy)/16.0)/log(1/zfactor);
	if (N < nscales) nscales = N;
    
	Dual_TVL1_optic_flow_multiscale(im1,im2, flv, flh, dimy, dimx, tau, lambda, theta, nscales, zfactor, nwarps, epsilon, 0);

    return;
}