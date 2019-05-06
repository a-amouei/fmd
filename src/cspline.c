/*
  cspline.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "cspline.h"

void spline_prepare_adv(double x[], double y[], int n, double yDD[])
/*
Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function, i.e., y_i = f ( x_i ), with
x_0 < x_1 < ... < x_(n-1), this routine returns an array yDD[0..n-1] that contains
the second derivatives of the interpolating function at the tabulated points x_i. The
boundary conditions are set for a natural spline, with zero second derivative on both boundaries.
*/
{
    int i,k;
    double p,sig,*u;

    u=malloc((n-1)*sizeof(double));
    yDD[0]=u[0]=0.0;
    for (i=1;i<=n-2;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*yDD[i-1]+2.0;
        yDD[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    yDD[n-1]=0.0;
    for (k=n-2;k>=0;k--)
        yDD[k]=yDD[k]*yDD[k+1]+u[k];
    free(u);
}

void spline_prepare(double h, double y[], int n, double yDD[])
/*
Given the array y[0..n-1] containing a tabulated function, i.e., y_i = f ( x_i ), with
x_i = h * i , this routine returns an array yDD[0..n-1] that contains
the second derivatives of the interpolating function at the points x_i. The
boundary conditions are set for a natural spline, with zero second derivative on both boundaries.
*/
{
    int i,k;
    double p,sig=0.5,*u;

    u=malloc((n-1)*sizeof(double));
    yDD[0]=u[0]=0.0;
    for (i=1;i<=n-2;i++)
    {
        p=sig*yDD[i-1]+2.0;
        yDD[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-2*y[i]+y[i-1])/h;
        u[i]=(3.0*u[i]/h-sig*u[i-1])/p;
    }
    yDD[n-1]=0.0;
    for (k=n-2;k>=0;k--)
        yDD[k]=yDD[k]*yDD[k+1]+u[k];
    free(u);
}

double spline_val(double h, double ya[], double yDDa[], double x)
/*  Given the array ya[0..n-1] which tabulates a function,
and given the array yDDa[0..n-1], which is the output from spline_prepare(), and given a value of
x, this routine returns a cubic-spline interpolated value. */
{
    int lo,hi;
    double b,a;

    lo=(int)(x / h);
    hi=lo+1;

    a=hi-x/h;
    b=1-a;
    return SPLINE_VAL(a,b,ya,lo,hi,yDDa,h);
}

double spline_deriv(double h, double ya[], double yDDa[], double x)
{
    int lo,hi;
    double b,a;

    lo=(int)(x / h);
    hi=lo+1;

    a=hi-x/h;
    b=1-a;
    return SPLINE_DERIV(a,b,ya,lo,hi,yDDa,h);
}

void spline_val_deriv(double h, double ya[], double yDDa[], double x, double *y, double *yD)
{
    int lo,hi;
    double b,a;

    lo=(int)(x / h);
    hi=lo+1;

    a=hi-x/h;
    b=1-a;
    *y=SPLINE_VAL(a,b,ya,lo,hi,yDDa,h);
    *yD=SPLINE_DERIV(a,b,ya,lo,hi,yDDa,h);
}
