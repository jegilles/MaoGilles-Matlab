
// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef DUAL_TVL1_OPTIC_FLOW_H
#define DUAL_TVL1_OPTIC_FLOW_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "mask.c"
#include "bicubic_interpolation.c"
#include "zoom.c"

#define MAX_ITERATIONS 300
#define PRESMOOTHING_SIGMA 0.8
#define GRAD_IS_ZERO 1E-10

/**
 * Implementation of the Zach, Pock and Bischof dual TV-L1 optic flow method
 *
 * see reference:
 *  [1] C. Zach, T. Pock and H. Bischof, "A Duality Based Approach for Realtime
 *      TV-L1 Optical Flow", In Proceedings of Pattern Recognition (DAGM),
 *      Heidelberg, Germany, pp. 214-223, 2007
 *
 *
 * Details on the total variation minimization scheme can be found in:
 *  [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 *      Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 **/


/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void Dual_TVL1_optic_flow(
		double *I0,           // source image
		double *I1,           // target image
		double *u1,           // x component of the optical flow
		double *u2,           // y component of the optical flow
		const int   nx,      // image width
		const int   ny,      // image height
		const double tau,     // time step
		const double lambda,  // weight parameter for the data term
		const double theta,   // weight parameter for (u - v)²
		const int   warps,   // number of warpings per scale
		const double epsilon, // tolerance for numerical convergence
		const bool  verbose  // enable/disable the verbose mode
		)
{
	const int   size = nx * ny;
	const double l_t = lambda * theta;

	size_t sf = sizeof(double);
	double *I1x    = malloc(size*sf);
	double *I1y    = xmalloc(size*sf);
	double *I1w    = xmalloc(size*sf);
	double *I1wx   = xmalloc(size*sf);
	double *I1wy   = xmalloc(size*sf);
	double *rho_c  = xmalloc(size*sf);
	double *v1     = xmalloc(size*sf);
	double *v2     = xmalloc(size*sf);
	double *p11    = xmalloc(size*sf);
	double *p12    = xmalloc(size*sf);
	double *p21    = xmalloc(size*sf);
	double *p22    = xmalloc(size*sf);
	double *div    = xmalloc(size*sf);
	double *grad   = xmalloc(size*sf);
	double *div_p1 = xmalloc(size*sf);
	double *div_p2 = xmalloc(size*sf);
	double *u1x    = xmalloc(size*sf);
	double *u1y    = xmalloc(size*sf);
	double *u2x    = xmalloc(size*sf);
	double *u2y    = xmalloc(size*sf);
	int i,warpings;

	centered_gradient(I1, I1x, I1y, nx, ny);

	// initialization of p
	for (i = 0; i < size; i++)
	{
		p11[i] = p12[i] = 0.0;
		p21[i] = p22[i] = 0.0;
	}

	for (warpings = 0; warpings < warps; warpings++)
	{
		// compute the warping of the target image and its derivatives
		bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
		bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
		bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

//#pragma omp parallel for
		for (i = 0; i < size; i++)
		{
			const double Ix2 = I1wx[i] * I1wx[i];
			const double Iy2 = I1wy[i] * I1wy[i];

			// store the |Grad(I1)|^2
			grad[i] = (Ix2 + Iy2);

			// compute the constant part of the rho function
			rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
						- I1wy[i] * u2[i] - I0[i]);
		}

		int n = 0;
		double error = INFINITY;
		while (error > epsilon * epsilon && n < MAX_ITERATIONS)
		{
			n++;
			// estimate the values of the variable (v1, v2)
			// (thresholding opterator TH)
//#pragma omp parallel for
			for (i = 0; i < size; i++)
			{
				const double rho = rho_c[i]
					+ (I1wx[i] * u1[i] + I1wy[i] * u2[i]);

				double d1, d2;

				if (rho < - l_t * grad[i])
				{
					d1 = l_t * I1wx[i];
					d2 = l_t * I1wy[i];
				}
				else
				{
					if (rho > l_t * grad[i])
					{
						d1 = -l_t * I1wx[i];
						d2 = -l_t * I1wy[i];
					}
					else
					{
						if (grad[i] < GRAD_IS_ZERO)
							d1 = d2 = 0;
						else
						{
							double fi = -rho/grad[i];
							d1 = fi * I1wx[i];
							d2 = fi * I1wy[i];
						}
					}
				}

				v1[i] = u1[i] + d1;
				v2[i] = u2[i] + d2;
			}

			// compute the divergence of the dual variable (p1, p2)
			divergence(p11, p12, div_p1, nx ,ny);
			divergence(p21, p22, div_p2, nx ,ny);

			// estimate the values of the optical flow (u1, u2)
			error = 0.0;
//#pragma omp parallel for reduction(+:error)
			for (i = 0; i < size; i++)
			{
				const double u1k = u1[i];
				const double u2k = u2[i];

				u1[i] = v1[i] + theta * div_p1[i];
				u2[i] = v2[i] + theta * div_p2[i];

				error += (u1[i] - u1k) * (u1[i] - u1k) +
					(u2[i] - u2k) * (u2[i] - u2k);
			}
			error /= size;

			// compute the gradient of the optical flow (Du1, Du2)
			forward_gradient(u1, u1x, u1y, nx ,ny);
			forward_gradient(u2, u2x, u2y, nx ,ny);

			// estimate the values of the dual variable (p1, p2)
//#pragma omp parallel for
			for (i = 0; i < size; i++)
			{
				const double taut = tau / theta;
				const double g1   = hypot(u1x[i], u1y[i]);
				const double g2   = hypot(u2x[i], u2y[i]);
				const double ng1  = 1.0 + taut * g1;
				const double ng2  = 1.0 + taut * g2;

				p11[i] = (p11[i] + taut * u1x[i]) / ng1;
				p12[i] = (p12[i] + taut * u1y[i]) / ng1;
				p21[i] = (p21[i] + taut * u2x[i]) / ng2;
				p22[i] = (p22[i] + taut * u2y[i]) / ng2;
			}
		}

		if (verbose)
			fprintf(stderr, "Warping: %d, "
					"Iterations: %d, "
					"Error: %f\n", warpings, n, error);
	}

	// delete allocated memory
	free(I1x);
	free(I1y);
	free(I1w);
	free(I1wx);
	free(I1wy);
	free(rho_c);
	free(v1);
	free(v2);
	free(p11);
	free(p12);
	free(p21);
	free(p22);
	free(div);
	free(grad);
	free(div_p1);
	free(div_p2);
	free(u1x);
	free(u1y);
	free(u2x);
	free(u2y);
}

/**
 *
 * Compute the max and min of an array
 *
 **/
static void getminmax(
	double *min,     // output min
	double *max,     // output max
	const double *x, // input array
	int n           // array size
)
{
	int i;
	*min = *max = x[0];
	for (i = 1; i < n; i++) {
		if (x[i] < *min)
			*min = x[i];
		if (x[i] > *max)
			*max = x[i];
	}
}

/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization(
		const double *I0,  // input image0
		const double *I1,  // input image1
		double *I0n,       // normalized output image0
		double *I1n,       // normalized output image1
		int size          // size of the image
		)
{
	double max0, max1, min0, min1;
	int i;

	// obtain the max and min of each image
	getminmax(&min0, &max0, I0, size);
	getminmax(&min1, &max1, I1, size);

	// obtain the max and min of both images
	const double max = (max0 > max1)? max0 : max1;
	const double min = (min0 < min1)? min0 : min1;
	const double den = max - min;

	if (den > 0)
		// normalize both images
		for (i = 0; i < size; i++)
		{
			I0n[i] = 255.0 * (I0[i] - min) / den;
			I1n[i] = 255.0 * (I1[i] - min) / den;
		}

	else
		// copy the original images
		for (i = 0; i < size; i++)
		{
			I0n[i] = I0[i];
			I1n[i] = I1[i];
		}
}


/**
 *
 * Function to compute the optical flow using multiple scales
 *
 **/
void Dual_TVL1_optic_flow_multiscale(
		double *I0,           // source image
		double *I1,           // target image
		double *u1,           // x component of the optical flow
		double *u2,           // y component of the optical flow
		const int   nxx,     // image width
		const int   nyy,     // image height
		const double tau,     // time step
		const double lambda,  // weight parameter for the data term
		const double theta,   // weight parameter for (u - v)²
		const int   nscales, // number of scales
		const double zfactor, // factor for building the image piramid
		const int   warps,   // number of warpings per scale
		const double epsilon, // tolerance for numerical convergence
		const bool  verbose  // enable/disable the verbose mode
)
{
	int size = nxx * nyy;

	// allocate memory for the pyramid structure
	double **I0s = xmalloc(nscales * sizeof(double*));
	double **I1s = xmalloc(nscales * sizeof(double*));
	double **u1s = xmalloc(nscales * sizeof(double*));
	double **u2s = xmalloc(nscales * sizeof(double*));
	int    *nx  = xmalloc(nscales * sizeof(int));
	int    *ny  = xmalloc(nscales * sizeof(int));
	int s,i;

	I0s[0] = xmalloc(size*sizeof(double));
	I1s[0] = xmalloc(size*sizeof(double));

	u1s[0] = u1;
	u2s[0] = u2;
	nx [0] = nxx;
	ny [0] = nyy;

	// normalize the images between 0 and 255
	image_normalization(I0, I1, I0s[0], I1s[0], size);

	// pre-smooth the original images
	gaussian(I0s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);
	gaussian(I1s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);

	// create the scales
	for (s = 1; s < nscales; s++)
	{
		zoom_size(nx[s-1], ny[s-1], &nx[s], &ny[s], zfactor);
		const int sizes = nx[s] * ny[s];

		// allocate memory
		I0s[s] = xmalloc(sizes*sizeof(double));
		I1s[s] = xmalloc(sizes*sizeof(double));
		u1s[s] = xmalloc(sizes*sizeof(double));
		u2s[s] = xmalloc(sizes*sizeof(double));

		// zoom in the images to create the pyramidal structure
		zoom_out(I0s[s-1], I0s[s], nx[s-1], ny[s-1], zfactor);
		zoom_out(I1s[s-1], I1s[s], nx[s-1], ny[s-1], zfactor);
	}

	// initialize the flow at the coarsest scale
	for (i = 0; i < nx[nscales-1] * ny[nscales-1]; i++)
		u1s[nscales-1][i] = u2s[nscales-1][i] = 0.0;

	// pyramidal structure for computing the optical flow
	for (s = nscales-1; s >= 0; s--)
	{
		if (verbose)
			fprintf(stderr, "Scale %d: %dx%d\n", s, nx[s], ny[s]);

		// compute the optical flow at the current scale
		Dual_TVL1_optic_flow(
				I0s[s], I1s[s], u1s[s], u2s[s], nx[s], ny[s],
				tau, lambda, theta, warps, epsilon, verbose
		);

		// if this was the last scale, finish now
		if (!s) break;

		// otherwise, upsample the optical flow

		// zoom the optical flow for the next finer scale
		zoom_in(u1s[s], u1s[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
		zoom_in(u2s[s], u2s[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);

		// scale the optical flow with the appropriate zoom factor
		for (i = 0; i < nx[s-1] * ny[s-1]; i++)
		{
			u1s[s-1][i] *= (double) 1.0 / zfactor;
			u2s[s-1][i] *= (double) 1.0 / zfactor;
		}
	}

	// delete allocated memory
	for (i = 1; i < nscales; i++)
	{
		free(I0s[i]);
		free(I1s[i]);
		free(u1s[i]);
		free(u2s[i]);
	}
	free(I0s[0]);
	free(I1s[0]);

	free(I0s);
	free(I1s);
	free(u1s);
	free(u2s);
	free(nx);
	free(ny);
}


#endif//DUAL_TVL1_OPTIC_FLOW_H
