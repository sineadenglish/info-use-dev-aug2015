#ifndef GROWTHFUN_H
#define GROWTHFUN_H

// define constants to be used in growth model (minimum and maximum size, number of levels for size, belief (P) and alpha)
const double	SIZEMIN = 1.0;
const double	SIZEMAX = 80.0; 
const int			SIZESTEPS = 80; 
const int			PSTEPS = 501; 
const int			ALPHASTEPS = 501; 

// define functions to be used in growth model (defined in growthfun.cpp file)
double maxval(double a[], int size);
double minval(double a[], int size);
double absol(double val);
double minval_2vec(double A[], double B[]);
double interpolate(double f[SIZESTEPS][PSTEPS], double x, double y);
double interpolate2(double f[PSTEPS][ALPHASTEPS], double x, double y);
double nearestval(double f[SIZESTEPS][PSTEPS], double x, double y);

#endif