// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2012, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>
// All rights reserved.


// 
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// this function is like "malloc", but it returns always a valid pointer
static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new)
		exit(fprintf(stderr, "out of memory\n"));
	return new;
}

// the type of the "getpixel" function
typedef double (*extension_operator_double)(double*, int, int, int, int);

// getpixel, with neumann boundary conditions
static double extend_double_image_constant(double *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

// compute the gradient and temporal derivative of the input image pair
static void compute_input_derivatives(double *Ex, double *Ey, double *Et,
		double *a, double *b, int w, int h)
{
    int i,j;
	extension_operator_double p = extend_double_image_constant;
	for (j = 0; j < h; j++)
	for (i = 0; i < w; i++) {
		Ey[j*w+i] = (1.0/4) * ( p(a,w,h, i, j+1) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i+1,j));
		Ex[j*w+i] = (1.0/4) * ( p(a,w,h, i+1, j) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i,j+1));
		Et[j*w+i] = (1.0/4) * ( p(b,w,h, i, j) - p(a,w,h, i,j)
				+ p(b,w,h, i+1, j) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j+1) - p(a,w,h, i+1,j+1));
	}
}

// compute a local average of a function "u"
static void compute_bar(double *ubar, double *u, int w, int h)
{
    int i,j;
	extension_operator_double p = extend_double_image_constant;
	for (j = 0; j < h; j++)
	for (i = 0; i < w; i++)
		ubar[j*w+i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));
}

// compute a sigle iteration of the classical Horn-Schunck method
static void hs_iteration(double *u, double *v,
		double *Ex, double *Ey, double *Et, int w, int h, double alpha)
{
    int i;
	double *ubar = xmalloc(w * h * sizeof(double));
	double *vbar = xmalloc(w * h * sizeof(double));
	compute_bar(ubar, u, w, h);
	compute_bar(vbar, v, w, h);
	for (i = 0; i < w*h; i++) {
		double t = Ex[i]*ubar[i] + Ey[i]*vbar[i] + Et[i];
		t /= alpha*alpha + Ex[i]*Ex[i] + Ey[i]*Ey[i];
		u[i] = ubar[i] - Ex[i] * t;
		v[i] = vbar[i] - Ey[i] * t;
	}
	free(ubar);
	free(vbar);
}

// run n iterations of the classical Horn-Schunck method
void hs(double *u, double *v, double *a, double *b, int w, int h,
		int n, double alpha)
{
	double *gx = xmalloc(w * h * sizeof(double));
	double *gy = xmalloc(w * h * sizeof(double));
	double *gt = xmalloc(w * h * sizeof(double));
    int i;
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	for (i = 0; i < w*h; i++)
		u[i] = v[i] = 0;
	for (i = 0; i < n; i++)
		hs_iteration(u, v, gx, gy, gt, w, h, alpha);
	free(gx);
	free(gy);
       	free(gt);
}