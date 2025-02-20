// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef BICUBIC_INTERPOLATION_C
#define BICUBIC_INTERPOLATION_C

#include <stdbool.h>

#define BOUNDARY_CONDITION 0
//0 Neumann
//1 Periodic
//2 Symmetric

/**
  *
  * Neumann boundary condition test
  *
**/
static int neumann_bc(int x, int nx, bool *out)
{
	if(x < 0)
	{
	    x = 0;
	    *out = true;
	}
	else if (x >= nx)
	{
	    x = nx - 1;
	    *out = true;
	}

	return x;
}

/**
  *
  * Periodic boundary condition test
  *
**/
static int periodic_bc(int x, int nx, bool *out)
{
	if(x < 0)
	{
		const int n   = 1 - (int)(x/(nx+1));
		const int ixx = x + n * nx;

		x =   ixx% nx;
		*out = true;
	}
	else if(x >= nx)
	{
		x = x % nx;
		*out = true;
	}

	return x;
}


/**
  *
  * Symmetric boundary condition test
  *
**/
static int symmetric_bc(int x, int nx, bool *out)
{
	if(x < 0)
	{
		const int borde = nx - 1;
		const int xx = -x;
		const int n  = (int)(xx/borde) % 2;

		if ( n ) x = borde - ( xx % borde );
		else x = xx % borde;
		*out = true;
	}
	else if ( x >= nx )
	{
		const int borde = nx - 1;
		const int n = (int)(x/borde) % 2;

		if ( n ) x = borde - ( x % borde );
		else x = x % borde;
		*out = true;
	}

	return x;
}


/**
  *
  * Cubic interpolation in one dimension
  *
**/
static double cubic_interpolation_cell (
	double v[4],  //interpolation points
	double x      //point to be interpolated
)
{
	return  v[1] + 0.5 * x * (v[2] - v[0] +
		x * (2.0 *  v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] +
		x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}


/**
  *
  * Bicubic interpolation in two dimensions
  *
**/
static double bicubic_interpolation_cell (
	double p[4][4], //array containing the interpolation points
	double x,       //x position to be interpolated
	double y        //y position to be interpolated
)
{
	double v[4];
	v[0] = cubic_interpolation_cell(p[0], y);
	v[1] = cubic_interpolation_cell(p[1], y);
	v[2] = cubic_interpolation_cell(p[2], y);
	v[3] = cubic_interpolation_cell(p[3], y);
	return cubic_interpolation_cell(v, x);
}

/**
  *
  * Compute the bicubic interpolation of a point in an image.
  * Detect if the point goes outside the image domain.
  *
**/
double bicubic_interpolation_at(
	const double *input, //image to be interpolated
	const double  uu,    //x component of the vector field
	const double  vv,    //y component of the vector field
	const int    nx,    //image width
	const int    ny,    //image height
	bool         border_out //if true, return zero outside the region
)
{
	const int sx = (uu < 0)? -1: 1;
	const int sy = (vv < 0)? -1: 1;

	int x, y, mx, my, dx, dy, ddx, ddy;
	bool out[1] = {false};

	//apply the corresponding boundary conditions
	switch(BOUNDARY_CONDITION) {

		case 0: x   = neumann_bc((int) uu, nx, out);
			y   = neumann_bc((int) vv, ny, out);
			mx  = neumann_bc((int) uu - sx, nx, out);
			my  = neumann_bc((int) vv - sx, ny, out);
			dx  = neumann_bc((int) uu + sx, nx, out);
			dy  = neumann_bc((int) vv + sy, ny, out);
			ddx = neumann_bc((int) uu + 2*sx, nx, out);
			ddy = neumann_bc((int) vv + 2*sy, ny, out);
			break;

		case 1: x   = periodic_bc((int) uu, nx, out);
			y   = periodic_bc((int) vv, ny, out);
			mx  = periodic_bc((int) uu - sx, nx, out);
			my  = periodic_bc((int) vv - sx, ny, out);
			dx  = periodic_bc((int) uu + sx, nx, out);
			dy  = periodic_bc((int) vv + sy, ny, out);
			ddx = periodic_bc((int) uu + 2*sx, nx, out);
			ddy = periodic_bc((int) vv + 2*sy, ny, out);
			break;

		case 2: x   = symmetric_bc((int) uu, nx, out);
			y   = symmetric_bc((int) vv, ny, out);
			mx  = symmetric_bc((int) uu - sx, nx, out);
			my  = symmetric_bc((int) vv - sx, ny, out);
			dx  = symmetric_bc((int) uu + sx, nx, out);
			dy  = symmetric_bc((int) vv + sy, ny, out);
			ddx = symmetric_bc((int) uu + 2*sx, nx, out);
			ddy = symmetric_bc((int) vv + 2*sy, ny, out);
			break;

		default:x   = neumann_bc((int) uu, nx, out);
			y   = neumann_bc((int) vv, ny, out);
			mx  = neumann_bc((int) uu - sx, nx, out);
			my  = neumann_bc((int) vv - sx, ny, out);
			dx  = neumann_bc((int) uu + sx, nx, out);
			dy  = neumann_bc((int) vv + sy, ny, out);
			ddx = neumann_bc((int) uu + 2*sx, nx, out);
			ddy = neumann_bc((int) vv + 2*sy, ny, out);
			break;
	}

	if(*out && border_out)
		return 0.0;

	else
	{
		//obtain the interpolation points of the image
		const double p11 = input[mx  + nx * my];
		const double p12 = input[x   + nx * my];
		const double p13 = input[dx  + nx * my];
		const double p14 = input[ddx + nx * my];

		const double p21 = input[mx  + nx * y];
		const double p22 = input[x   + nx * y];
		const double p23 = input[dx  + nx * y];
		const double p24 = input[ddx + nx * y];

		const double p31 = input[mx  + nx * dy];
		const double p32 = input[x   + nx * dy];
		const double p33 = input[dx  + nx * dy];
		const double p34 = input[ddx + nx * dy];

		const double p41 = input[mx  + nx * ddy];
		const double p42 = input[x   + nx * ddy];
		const double p43 = input[dx  + nx * ddy];
		const double p44 = input[ddx + nx * ddy];

		//create array
		double pol[4][4] = {
			{p11, p21, p31, p41},
			{p12, p22, p32, p42},
			{p13, p23, p33, p43},
			{p14, p24, p34, p44}
		};

		//return interpolation
		return bicubic_interpolation_cell(pol, uu-x, vv-y);
	}
}


/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
void bicubic_interpolation_warp(
	const double *input,     // image to be warped
	const double *u,         // x component of the vector field
	const double *v,         // y component of the vector field
	double       *output,    // image warped with bicubic interpolation
	const int    nx,        // image width
	const int    ny,        // image height
	bool         border_out // if true, put zeros outside the region
)
{
	int i,j;
//	#pragma omp parallel for
	for(i = 0; i < ny; i++)
		for(j = 0; j < nx; j++)
		{
			const int   p  = i * nx + j;
			const double uu = (double) (j + u[p]);
			const double vv = (double) (i + v[p]);

			// obtain the bicubic interpolation at position (uu, vv)
			output[p] = bicubic_interpolation_at(input,
					uu, vv, nx, ny, border_out);
		}
}


#endif//BICUBIC_INTERPOLATION_C
