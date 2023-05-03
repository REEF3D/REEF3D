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

#include"particle_pls.h"
#include"fdm.h"
#include<math.h>

double particle_pls::phipol(lexer *p,fdm* a, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
		if(xp/dx-0.5<0.0)
		--i;
    j=int((yp)/dx-0.5);
		if(yp/dx-0.5<0.0)
		--j;
    k=int((zp)/dx-0.5);
		if(zp/dx-0.5<0.0)
		--k;

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    if(ipolval==1)
    value =  lint(a->phi,i,j,k,wa,wb,wc);

    if(ipolval==2)
    value =  tricubic(p,a,a->phi,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double particle_pls::upol(lexer *p,fdm* a, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-1.0);
		if(xp/dx-1.0<0.0)
		--i;
		
    j=int((yp)/dx-0.5);
		if(yp/dx-0.5<0.0)
		--j;
		
    k=int((zp)/dx-0.5);
		if(zp/dx-0.5<0.0)
		--k;

    wa=((double(i) + 2.0)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    if(ipolval==1)
    value =  lint(a->u,i,j,k,wa,wb,wc);

    if(ipolval==2)
    value =  tricubic(p,a,a->u,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double particle_pls::vpol(lexer *p,fdm* a, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
		if(xp/dx-0.5<0.0)
		--i;
		
    j=int((yp)/dx-1.0);
		if(yp/dx-1.0<0.0)
		--j;
		
    k=int((zp)/dx-0.5);
		if(zp/dx-0.5<0.0)
		--k;

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 2.0)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    if(ipolval==1)
    value =  lint(a->v,i,j,k,wa,wb,wc);

    if(ipolval==2)
    value =  tricubic(p,a,a->v,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double particle_pls::wpol(lexer *p,fdm* a, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
		if(xp/dx-0.5<0.0)
		--i;
		
    j=int((yp)/dx-0.5);
		if(yp/dx-0.5<0.0)
		--j;
		
    k=int((zp)/dx-1.0);
		if(zp/dx-1.0<0.0)
		--k;

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 2.0)-zp/dx);

    if(ipolval==1)
    value =  lint(a->w,i,j,k,wa,wb,wc);

    if(ipolval==2)
    value =  tricubic(p,a,a->w,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double particle_pls::lint(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
pip=4;

    x1 = wa*f(i,j,k)   + (1.0-wa)*f(i+1,j,k);
    x2 = wa*f(i,j+1,k) + (1.0-wa)*f(i+1,j+1,k);

    x3 = wa*f(i,j,k+1)   + (1.0-wa)*f(i+1,j,k+1);
    x4 = wa*f(i,j+1,k+1) + (1.0-wa)*f(i+1,j+1,k+1);

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;

pip=0;
 return value;

}

double particle_pls::cint(double wx, double f0, double f1, double f2, double f3)
{
    di0 = (f2-f0)*0.5;
    di1 = (f3-f1)*0.5;
    df=(f2-f1);

    if(fabs(df)<1.0e-15)
    di0=di1=0.0;

    if(di0*df<0.0)
    di0=0.0;

    if(di1*df<0.0)
    di1=0.0;

    a0 = f1;
    a1 = di0;
    a2 = 3.0*df - 2.0*di0 - di1;
    a3 = di0 + di1 - 2.0*df;

    x1 = a0 + a1*(1.0-wx) + a2*pow(1.0-wx,2.0) + a3*pow(1.0-wx,3.0);


return x1;
}

double particle_pls::tricubic(lexer *p,fdm* a,field& f, int& i,int& j, int& k, double wx, double wy, double wz)
{
pip=4;
    i0=j0=k0=1;
    i3=j3=k3=2;

    if(minboundcheck(p,i,j,k,1)==0)
    i0=j0=k0=0;

    if(maxboundcheck(p,i,j,k,1)==0)
    i3=j3=k3=0;

	value =  cint(wz, cint(wx, cint(wy, f(i-i0,j-j0,k-k0),f(i-i0,j,k-k0),
                              f(i-i0,j+1,k-k0),f(i-i0,j+j3,k-k0)),
                              cint(wy, f(i,j-j0,k-k0),f(i,j,k-k0),
                              f(i,j+1,k-k0),f(i,j+j3,k-k0)),
                              cint(wy, f(i+1,j-j0,k-k0),f(i+1,j,k-k0),
                              f(i+1,j+1,k-k0),f(i+1,j+j3,k-k0)),
                              cint(wy, f(i+i3,j-j0,k-k0),f(i+i3,j,k-k0),
                              f(i+i3,j+1,k-k0),f(i+i3,j+j3,k-k0)) ),
                     cint(wx, cint(wy, f(i-i0,j-j0,k),f(i-i0,j,k),
                              f(i-i0,j+1,k),f(i-i0,j+j3,k)),
                              cint(wy, f(i,j-j0,k),f(i,j,k),
                              f(i,j+1,k),f(i,j+j3,k)),
                              cint(wy, f(i+1,j-j0,k),f(i+1,j,k),
                              f(i+1,j+1,k),f(i+1,j+j3,k)),
                              cint(wy, f(i+i3,j-j0,k),f(i+i3,j,k),
                              f(i+i3,j+1,k),f(i+i3,j+j3,k)) ),
                     cint(wx, cint(wy, f(i-i0,j-j0,k+1),f(i-i0,j,k+1),
                              f(i-i0,j+1,k+1),f(i-i0,j+j3,k+1)),
                              cint(wy, f(i,j-j0,k+1),f(i,j,k+1),
                              f(i,j+1,k+1),f(i,j+j3,k+1)),
                              cint(wy, f(i+1,j-j0,k+1),f(i+1,j,k+1),
                              f(i+1,j+1,k+1),f(i+1,j+j3,k+1)),
                              cint(wy, f(i+i3,j-j0,k+1),f(i+i3,j,k+1),
                              f(i+i3,j+1,k+1),f(i+i3,j+j3,k+1)) ),
                     cint(wx, cint(wy, f(i-i0,j-j0,k+k3),f(i-i0,j,k+k3),
                              f(i-i0,j+1,k+k3),f(i-i0,j+j3,k+k3)),
                              cint(wy, f(i,j-j0,k+k3),f(i,j,k+k3),
                              f(i,j+1,k+k3),f(i,j+j3,k+k3)),
                              cint(wy, f(i+1,j-j0,k+k3),f(i+1,j,k+k3),
                              f(i+1,j+1,k+k3),f(i+1,j+j3,k+k3)),
                              cint(wy, f(i+i3,j-j0,k+k3),f(i+i3,j,k+k3),
                              f(i+i3,j+1,k+k3),f(i+i3,j+j3,k+k3)) ) );
pip=0;

    return value;
}


