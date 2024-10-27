// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef MASK_C
#define MASK_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "xmalloc.c"


#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING


/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/


/**
 *
 * Function to compute the divergence with backward differences
 * (see [2] for details)
 *
 **/
void divergence(
		const double *v1, // x component of the vector field
		const double *v2, // y component of the vector field
		double *div,      // output divergence
		const int nx,    // image width
		const int ny     // image height
	       )
{
		int i,j;
	// compute the divergence on the central body of the image
//#pragma omp parallel for schedule(dynamic)
	for (i = 1; i < ny-1; i++)
	{
		for(j = 1; j < nx-1; j++)
		{
			const int p  = i * nx + j;
			const int p1 = p - 1;
			const int p2 = p - nx;

			const double v1x = v1[p] - v1[p1];
			const double v2y = v2[p] - v2[p2];

			div[p] = v1x + v2y;
		}
	}

	// compute the divergence on the first and last rows
	for (j = 1; j < nx-1; j++)
	{
		const int p = (ny-1) * nx + j;

		div[j] = v1[j] - v1[j-1] + v2[j];
		div[p] = v1[p] - v1[p-1] - v2[p-nx];
	}

	// compute the divergence on the first and last columns
	for (i = 1; i < ny-1; i++)
	{
		const int p1 = i * nx;
		const int p2 = (i+1) * nx - 1;

		div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
		div[p2] = -v1[p2-1] + v2[p2] - v2[p2 - nx];

	}

	div[0]         =  v1[0] + v2[0];
	div[nx-1]      = -v1[nx - 2] + v2[nx - 1];
	div[(ny-1)*nx] =  v1[(ny-1)*nx] - v2[(ny-2)*nx];
	div[ny*nx-1]   = -v1[ny*nx - 2] - v2[(ny-1)*nx - 1];
}


/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void forward_gradient(
		const double *f, //input image
		double *fx,      //computed x derivative
		double *fy,      //computed y derivative
		const int nx,   //image width
		const int ny    //image height
		)
{
	int i,j;
	// compute the gradient on the central body of the image
//#pragma omp parallel for schedule(dynamic)
	for (i = 0; i < ny-1; i++)
	{
		for(j = 0; j < nx-1; j++)
		{
			const int p  = i * nx + j;
			const int p1 = p + 1;
			const int p2 = p + nx;

			fx[p] = f[p1] - f[p];
			fy[p] = f[p2] - f[p];
		}
	}

	// compute the gradient on the last row
	for (j = 0; j < nx-1; j++)
	{
		const int p = (ny-1) * nx + j;

		fx[p] = f[p+1] - f[p];
		fy[p] = 0;
	}

	// compute the gradient on the last column
	for (i = 1; i < ny; i++)
	{
		const int p = i * nx-1;

		fx[p] = 0;
		fy[p] = f[p+nx] - f[p];
	}

	fx[ny * nx - 1] = 0;
	fy[ny * nx - 1] = 0;
}


/**
 *
 * Function to compute the gradient with centered differences
 *
 **/
void centered_gradient(
		const double *input,  //input image
		double *dx,           //computed x derivative
		double *dy,           //computed y derivative
		const int nx,        //image width
		const int ny         //image height
		)
{
	int i,j;
	// compute the gradient on the center body of the image
//#pragma omp parallel for schedule(dynamic)
	for (i = 1; i < ny-1; i++)
	{
		for(j = 1; j < nx-1; j++)
		{
			const int k = i * nx + j;
			dx[k] = 0.5*(input[k+1] - input[k-1]);
			dy[k] = 0.5*(input[k+nx] - input[k-nx]);
		}
	}

	// compute the gradient on the first and last rows
	for (j = 1; j < nx-1; j++)
	{
		dx[j] = 0.5*(input[j+1] - input[j-1]);
		dy[j] = 0.5*(input[j+nx] - input[j]);

		const int k = (ny - 1) * nx + j;

		dx[k] = 0.5*(input[k+1] - input[k-1]);
		dy[k] = 0.5*(input[k] - input[k-nx]);
	}

	// compute the gradient on the first and last columns
	for(i = 1; i < ny-1; i++)
	{
		const int p = i * nx;
		dx[p] = 0.5*(input[p+1] - input[p]);
		dy[p] = 0.5*(input[p+nx] - input[p-nx]);

		const int k = (i+1) * nx - 1;

		dx[k] = 0.5*(input[k] - input[k-1]);
		dy[k] = 0.5*(input[k+nx] - input[k-nx]);
	}

	// compute the gradient at the four corners
	dx[0] = 0.5*(input[1] - input[0]);
	dy[0] = 0.5*(input[nx] - input[0]);

	dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
	dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

	dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
	dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

	dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
	dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}


/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
void gaussian(
	double *I,             // input/output image
	const int xdim,       // image width
	const int ydim,       // image height
	const double sigma    // Gaussian sigma
)
{
	const int boundary_condition = DEFAULT_BOUNDARY_CONDITION;
	const int window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

	const double den  = 2*sigma*sigma;
	const int   size = (int) (window_size * sigma) + 1 ;
	const int   bdx  = xdim + size;
	const int   bdy  = ydim + size;
	int i,j,k;

	if (boundary_condition && size > xdim) {
		fprintf(stderr, "GaussianSmooth: sigma too large\n");
		abort();
	}

	// compute the coefficients of the 1D convolution kernel
	double B[size];
	for(i = 0; i < size; i++)
		B[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);

	// normalize the 1D convolution kernel
	double norm = 0;
	for(i = 0; i < size; i++)
		norm += B[i];
	norm *= 2;
	norm -= B[0];
	for(i = 0; i < size; i++)
		B[i] /= norm;

	// convolution of each line of the input image
	double *R = xmalloc((size + xdim + size)*sizeof*R);

	for (k = 0; k < ydim; k++)
	{
		for (i = size; i < bdx; i++)
			R[i] = I[k * xdim + i - size];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for(i = 0, j = bdx; i < size; i++, j++)
				R[i] = R[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + size-i];
				R[j] = I[k * xdim + xdim-i-1];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + xdim-size+i];
				R[j] = I[k * xdim + i];
			}
			break;
		}

		for (i = size; i < bdx; i++)
		{
			double sum = B[0] * R[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * ( R[i-j] + R[i+j] );
			I[k * xdim + i - size] = sum;
		}
	}

	// convolution of each column of the input image
	double *T = xmalloc((size + ydim + size)*sizeof*T);

	for (k = 0; k < xdim; k++)
	{
		for (i = size; i < bdy; i++)
			T[i] = I[(i - size) * xdim + k];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for (i = 0, j = bdy; i < size; i++, j++)
				T[i] = T[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for (i = 0, j = bdy; i < size; i++, j++) {
				T[i] = I[(size-i) * xdim + k];
				T[j] = I[(ydim-i-1) * xdim + k];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for( i = 0, j = bdx; i < size; i++, j++) {
				T[i] = I[(ydim-size+i) * xdim + k];
				T[j] = I[i * xdim + k];
			}
			break;
		}

		for (i = size; i < bdy; i++)
		{
			double sum = B[0] * T[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * (T[i-j] + T[i+j]);
			I[(i - size) * xdim + k] = sum;
		}
	}

	free(R);
	free(T);
}


#endif//MASK_C
