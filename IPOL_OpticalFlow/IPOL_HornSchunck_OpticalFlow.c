/**************************************************************************
 *
 *  This code interface the Pyramidal Horn-Schunck OpticalFlow function 
 *  available at IPOL.
 *  The original C source code was modified to handle double instead of 
 *  float for an easier compatibility with Matlab
 *
 *************************************************************************/

#include <matrix.h>
#include <mex.h>  
#include "xmalloc.c"
#include "horn_schunck_pyramidal.c"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 //declare variables
    mxArray *in_im1, *in_im2, *out_flv, *out_flh,*in_alpha,*in_zfactor,*in_tol,*in_nwarps,*in_nscales,*in_maxiter;
    const mwSize *dims;   

    double *im1,*im2,*flv,*flh,alpha,zfactor,tol,N;
    int dimx,dimy,nwarps,nscales,maxiter;
    
 //associate Matlab inputs to duplicate C variables
    in_im1 = mxDuplicateArray(prhs[0]);
    in_im2 = mxDuplicateArray(prhs[1]);
    in_alpha = mxDuplicateArray(prhs[2]);
    in_zfactor = mxDuplicateArray(prhs[3]);
    in_tol = mxDuplicateArray(prhs[4]);
    in_nwarps = mxDuplicateArray(prhs[5]);
    in_nscales = mxDuplicateArray(prhs[6]);
    in_maxiter = mxDuplicateArray(prhs[6]);

//Get the dimensions of needed variables
    dims = mxGetDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    dims = mxGetDimensions(prhs[1]);
    if ((dimx != (int)dims[1])||(dimy != (int)dims[0])) mexErrMsgTxt("The two input images must be of the same size!");

//associate Matlab outputs to their corresponding C variables
    out_flv = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    out_flh = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//get memory pointers of certains duplicate variables
    im1 = mxGetPr(in_im1);
    im2 = mxGetPr(in_im2);
    flv = mxGetPr(out_flv); //vertical component of the flow
    flh = mxGetPr(out_flh); //horizontal component of the flow

//Get scalars of other duplicate variables
    alpha = mxGetScalar(in_alpha);
    zfactor = mxGetScalar(in_zfactor);
    tol = mxGetScalar(in_tol);
    nwarps = (int)mxGetScalar(in_nwarps);
    nscales = (int)mxGetScalar(in_nscales);
    maxiter = (int)mxGetScalar(in_maxiter);
    
 //run the TV-L1 algorithm originally coded in C
    N = 1+log(sqrt(dimx*dimx+dimy*dimy)/16.0)/log(1/zfactor);
	if (N < nscales) nscales = N;
    
	//we call the C function
    horn_schunck_pyramidal(im1, im2, flv, flh, dimy, dimx, alpha, nscales, zfactor, nwarps, tol, maxiter,0);

    return;
}
