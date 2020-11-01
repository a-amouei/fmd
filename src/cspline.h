/*
  cspline.h: This file is part of Free Molecular Dynamics

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

#ifndef CSPLINE_H
#define CSPLINE_H

#include "config.h"
#include "types.h"

#define SPLINE_VAL(a,b,y,lo,hi,yDD,h) ((a)*(y)[lo]+(b)*(y)[hi]+(((a)*(a)*(a)-(a))*(yDD)[lo]+((b)*(b)*(b)-(b))*(yDD)[hi])*((h)*(h))/6.0)
#define SPLINE_DERIV(a,b,y,lo,hi,yDD,h) (((y)[hi]-(y)[lo])/(h)+((3*(b)*(b)-1)*(yDD)[hi]-(3*(a)*(a)-1)*(yDD)[lo])*(h)/6.0)

#define FIND_klo_AND_ETC(list, x, n)                                   \
{                                                                      \
    klo=0;                                                             \
    khi=n-1;                                                           \
    while (khi-klo > 1)                                                \
    {                                                                  \
        k = (khi+klo) >> 1;                                            \
        if (list[k] > x)                                               \
            khi = k;                                                   \
        else                                                           \
            klo = k;                                                   \
    }                                                                  \
    h = list[khi] - list[klo];                                         \
    a=(list[khi]-x)/h;                                                 \
    b=(x-list[klo])/h;                                                 \
}

void spline_prepare_adv(fmd_real_t x[], fmd_real_t y[], int n, fmd_real_t yDD[]);
void spline_prepare(fmd_real_t h, fmd_real_t y[], int n, fmd_real_t yDD[]);
fmd_real_t spline_val(fmd_real_t h, fmd_real_t ya[], fmd_real_t yDDa[], fmd_real_t x);
fmd_real_t spline_deriv(fmd_real_t h, fmd_real_t ya[], fmd_real_t yDDa[], fmd_real_t x);
void spline_val_deriv(fmd_real_t h, fmd_real_t ya[], fmd_real_t yDDa[], fmd_real_t x, fmd_real_t *y, fmd_real_t *yD);

#endif /* CSPLINE_H */
