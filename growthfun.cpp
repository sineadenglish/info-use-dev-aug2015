// Define all constants and functions

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <vector> 
#include <iostream>
#include "growthfun.h"

using namespace std;

// return maximum value in an array or vector
double maxval(double a[], int size) 
{
	double max_val = a[0];    // initial value starts at the value of the first element instead of zero
	for (int i=1; i<size; i++) {
		if (a[i] > max_val) {
			max_val = a[i];
		}
	}
	return max_val;
}

// return minimum value in an array or vector
double minval(double a[], int size) 
{
	double min_val = a[0];    // initial value starts at the value of the first element instead of zero
	for (int i=1; i<size; i++) {
		if (a[i] < min_val) {
			min_val = a[i];
		}
	}
	return min_val;
}

// return the absolute value 
double absol(double val)
{
  if (val < 0.0)
	{return -val;}
  else        
	{return val;}
}     

// for two vectors (A and B), return the value in A at the location where B is minimum (ignoring NAN value of A)
double minval_2vec(double A[], double B[])
{
	// incorporate NA value into B:
	for(int i=0; i<4; i++) 
	{
		B[i] = B[i]*A[i]/A[i];
	}
	
	double min=B[0];
	int min_loc = 0;
	
	for (int i=1;i<4;i++)
	{
		if (B[i]<min | isnan(min)) 
		{
			min=B[i];
			min_loc=i;
		}    
	}
	
	return(A[min_loc]);
}


// have written two interpolate functions for the different sizes of 2D arrays (Note Andy: messy code but wasn't clever enough to do this independent of the constants)
// interpolate f(x,y): interpolate between two cells across a 2D array (for all size-belief states)
double interpolate(double f[SIZESTEPS][PSTEPS], double x, double y)
{
	int	ix,iy;						
	double xc,dx,yc,dy;	
	
	xc = (x - SIZEMIN) * (SIZESTEPS-1.0) / (SIZEMAX	- SIZEMIN);
	yc = (y - 0.0) * (PSTEPS-1.0) / (1 - 0.0);
	ix = floor(xc);
	iy = floor(yc);						
	dx = xc - ix;
	dy = yc - iy;
	
	if (ix == SIZESTEPS && iy == PSTEPS) return (f[ix][iy]);
	else if (ix < SIZESTEPS && iy == PSTEPS) return ((1-dx)*f[ix][iy] + dx*f[ix+1][iy]);
	else if (ix == SIZESTEPS && iy < PSTEPS) return ((1-dy)*f[ix][iy] + dy*f[ix][iy+1]);
	else return ((1-dx)*(1-dy)*f[ix][iy] + (1-dx)*dy*f[ix][iy+1] + dx*(1-dy)*f[ix+1][iy] + dx*dy*f[ix+1][iy+1]);	
}

// interpolate2 f(x,y): interpolate between two cells across a 2D array (for all belief-alpha states)
double interpolate2(double f[PSTEPS][ALPHASTEPS], double x, double y)
{
	int	ix,iy;						
	double xc,dx,yc,dy;	
	
	xc = (x - 0.0) * (PSTEPS-1.0) / (1 - 0.0);
	yc = (y - 0.0) * (ALPHASTEPS-1.0) / (1 - 0.0);
	ix = floor(xc);
	iy = floor(yc);						
	dx = xc - ix;
	dy = yc - iy;
	
	if (ix == PSTEPS && iy == ALPHASTEPS) return (f[ix][iy]);
	else if (ix < PSTEPS && iy == ALPHASTEPS) return ((1-dx)*f[ix][iy] + dx*f[ix+1][iy]);
	else if (ix == PSTEPS && iy < ALPHASTEPS) return ((1-dy)*f[ix][iy] + dy*f[ix][iy+1]);
	else return ((1-dx)*(1-dy)*f[ix][iy] + (1-dx)*dy*f[ix][iy+1] + dx*(1-dy)*f[ix+1][iy] + dx*dy*f[ix+1][iy+1]);	
}


// choose nearest non-NA value in grid
double nearestval(double f[SIZESTEPS][PSTEPS], double x, double y)
{
	int	ix,iy,incx,incy;	// integer values of x,y coordinates; and increments if looking in-between					
	double xc,dx,yc,dy;		// continuous values of x,y coordinates
	double bl,br,tl,tr;		// values of f[] at vertices of cell of interest (so can find nearest one)
	double allval[4], distvec[4]; // vector of vertix values and distance of value from each vertix
	
	xc = (x - SIZEMIN) * (SIZESTEPS-1.0) / (SIZEMAX	- SIZEMIN);
	yc = (y - 0.0) * (PSTEPS-1.0) / (1 - 0.0);
	ix = floor(xc);
	iy = floor(yc);						
	dx = xc - ix;
	dy = yc - iy;
	
	// get value at all four vertices
	incx = 1;
	incy = 1;
	if(xc==SIZEMAX) incx=0;
	if(yc==1.0) incy=0;
	
	bl = f[ix][iy];
	tl = f[ix][iy+incy];
	tr = f[ix+incx][iy+incy];
	br = f[ix+incx][iy];
	
	allval[0] = bl;
	allval[1] = tl;
	allval[2] = tr;
	allval[3] = br;

	double incx_d = incx; // recast as double for calculation below
	double incy_d = incy;
	
	distvec[0] = pow((pow(dx,2.0)+pow(dy,2.0)), 0.5);
	distvec[1] = pow((pow(dx,2.0)+pow((incy_d-dy),2.0)), 0.5);
	distvec[2] = pow((pow((incx_d-dx),2.0)+pow((incy_d-dy),2.0)), 0.5);
	distvec[3] = pow((pow((incx_d-dx),2.0)+pow(dy,2.0)), 0.5);
	
	// output the (non-NA) value at the vertix which is closest to the interpolated point
	return(minval_2vec(allval, distvec));
	
}



