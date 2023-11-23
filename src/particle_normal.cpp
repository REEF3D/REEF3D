/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_f::normal(fdm* a, double& xp, double& yp, double& zp, double& value)
{

    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
    j=int((yp)/dy-0.5);
    k=int((zp)/dz-0.5);


    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dy);
    wc=((double(k) + 1.5)-zp/dz);

    wx = wa*dx;
    wy = wb*dy;
    wz = wc*dz;

// i
    val1=wc*(wb*a->phi(i,j,k) + (1.0-wb)*a->phi(i,j+1,k))     + (1.0-wc)*(wb*a->phi(i,j,k+1) + (1.0-wb)*a->phi(i,j+1,k+1));
    val2=wc*(wb*a->phi(i+1,j,k) + (1.0-wb)*a->phi(i+1,j+1,k)) + (1.0-wc)*(wb*a->phi(i+1,j,k+1) + (1.0-wb)*a->phi(i+1,j+1,k+1));

    if(wa>0.8)
    di=(val2-value)/(wx);

    if(wa<0.2)
    di=(value-val1)/(dx-wx);

    if(wa>=0.2 && wa<=0.8)
    di=(val2*(dx-wx)*(dx-wx) - val1*wx*wx  + value*(wx*wx-(dx-wx)*(dx-wx)))/(wx*(dx-wx)*(wx+(dx-wx)));

// j
    val1=wc*(wa*a->phi(i,j,k) + (1.0-wa)*a->phi(i+1,j,k))     + (1.0-wc)*(wa*a->phi(i,j,k+1) + (1.0-wa)*a->phi(i+1,j,k+1));
    val2=wc*(wa*a->phi(i,j+1,k) + (1.0-wa)*a->phi(i+1,j+1,k)) + (1.0-wc)*(wa*a->phi(i,j+1,k+1) + (1.0-wa)*a->phi(i+1,j+1,k+1));


    if(wb>0.8)
    dj=(val2-value)/(wy);

    if(wb<0.2)
    dj=(value-val1)/((dy-wy));

    if(wb>=0.2 && wb<=0.8)
    dj=(val2*(dy-wy)*(dy-wy) -val1*wy*wy + value*(wy*wy-(dy-wy)*(dy-wy)))/(wy*(dy-wy)*(wy+(dy-wy)));

// k
    val1=wb*(wa*a->phi(i,j,k) + (1.0-wa)*a->phi(i+1,j,k))     + (1.0-wb)*(wa*a->phi(i,j+1,k) + (1.0-wa)*a->phi(i+1,j+1,k));
    val2=wb*(wa*a->phi(i,j,k+1) + (1.0-wa)*a->phi(i+1,j,k+1)) + (1.0-wb)*(wa*a->phi(i,j+1,k+1) + (1.0-wa)*a->phi(i+1,j+1,k+1));

    if(wc>0.8)
    dk=(val2-value)/(wz);

    if(wc<0.2)
    dk=(value-val1)/(dz-wz);

    if(wc>=0.2 && wc<=0.8)
    dk=(val2*(dz-wz)*(dz-wz) -val1*wz*wz + value*(wz*wz-(dz-wz)*(dz-wz)))/(wz*(dz-wz)*(wz+(dz-wz)));

    dnorm=sqrt(di*di + dj*dj + dk*dk);

    nvec[0]=di/(dnorm>1.0e-15?dnorm:1.0e20);
    nvec[1]=dj/(dnorm>1.0e-15?dnorm:1.0e20);
    nvec[2]=dk/(dnorm>1.0e-15?dnorm:1.0e20);

    i=ii;
    j=jj;
    k=kk;
}

void particle_f::normreg(fdm* a, int ii, int jj, int kk)
{
    i=ii;
    j=jj;
    k=kk;

    di = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(2.0*dx);
    dj = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(2.0*dy);
    dk = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(2.0*dz);

    dnorm=sqrt(di*di + dj*dj + dk*dk);

    nvec[0]=di/(dnorm>1.0e-15?dnorm:1.0e20);
    nvec[1]=dj/(dnorm>1.0e-15?dnorm:1.0e20);
    nvec[2]=dk/(dnorm>1.0e-15?dnorm:1.0e20);
}
