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

#include"interpolation.h"
#include"field.h"
#include"lexer.h"


double interpolation::tricubic4a(field& f, int& i,int& j, int& k, double wx, double wy, double wz)
{
    int i0,j0,k0,i3,j3,k3;
    
pip=4;
    i0=j0=k0=1;
    i3=j3=k3=2;

    if(minboundcheck(p,i,j,k,1)==0)
    i0=j0=k0=0;

    if(maxboundcheck(p,i,j,k,1)==0)
    i3=j3=k3=0;

	value =  cint4a(wz, cint4a(wx, cint4a(wy, f(i-i0,j-j0,k-k0),f(i-i0,j,k-k0),
                              f(i-i0,j+1,k-k0),f(i-i0,j+j3,k-k0)),
                              cint4a(wy, f(i,j-j0,k-k0),f(i,j,k-k0),
                              f(i,j+1,k-k0),f(i,j+j3,k-k0)),
                              cint4a(wy, f(i+1,j-j0,k-k0),f(i+1,j,k-k0),
                              f(i+1,j+1,k-k0),f(i+1,j+j3,k-k0)),
                              cint4a(wy, f(i+i3,j-j0,k-k0),f(i+i3,j,k-k0),
                              f(i+i3,j+1,k-k0),f(i+i3,j+j3,k-k0)) ),
                     cint4a(wx, cint4a(wy, f(i-i0,j-j0,k),f(i-i0,j,k),
                              f(i-i0,j+1,k),f(i-i0,j+j3,k)),
                              cint4a(wy, f(i,j-j0,k),f(i,j,k),
                              f(i,j+1,k),f(i,j+j3,k)),
                              cint4a(wy, f(i+1,j-j0,k),f(i+1,j,k),
                              f(i+1,j+1,k),f(i+1,j+j3,k)),
                              cint4a(wy, f(i+i3,j-j0,k),f(i+i3,j,k),
                              f(i+i3,j+1,k),f(i+i3,j+j3,k)) ),
                     cint4a(wx, cint4a(wy, f(i-i0,j-j0,k+1),f(i-i0,j,k+1),
                              f(i-i0,j+1,k+1),f(i-i0,j+j3,k+1)),
                              cint4a(wy, f(i,j-j0,k+1),f(i,j,k+1),
                              f(i,j+1,k+1),f(i,j+j3,k+1)),
                              cint4a(wy, f(i+1,j-j0,k+1),f(i+1,j,k+1),
                              f(i+1,j+1,k+1),f(i+1,j+j3,k+1)),
                              cint4a(wy, f(i+i3,j-j0,k+1),f(i+i3,j,k+1),
                              f(i+i3,j+1,k+1),f(i+i3,j+j3,k+1)) ),
                     cint4a(wx, cint4a(wy, f(i-i0,j-j0,k+k3),f(i-i0,j,k+k3),
                              f(i-i0,j+1,k+k3),f(i-i0,j+j3,k+k3)),
                              cint4a(wy, f(i,j-j0,k+k3),f(i,j,k+k3),
                              f(i,j+1,k+k3),f(i,j+j3,k+k3)),
                              cint4a(wy, f(i+1,j-j0,k+k3),f(i+1,j,k+k3),
                              f(i+1,j+1,k+k3),f(i+1,j+j3,k+k3)),
                              cint4a(wy, f(i+i3,j-j0,k+k3),f(i+i3,j,k+k3),
                              f(i+i3,j+1,k+k3),f(i+i3,j+j3,k+k3)) ) );
pip=0;

    return value;
}

double interpolation::cctripol4_a(fdm* a,field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int(xp/p->DXP[IP]-0.5);
		if(xp/p->DXP[IP]-0.5<0.0)
		--i;
    j=int(yp/p->DYP[JP]-0.5);
		if(yp/p->DYP[JP]-0.5<0.0)
		--j;
    k=int(zp/p->DZP[KP]-0.5);
		if(zp/p->DZP[KP]-0.5<0.0)
		--k;

    wa=p->XP[IP1]-xp/p->DXP[IP];
    wb=p->YP[JP2]-yp/p->DYP[JP];
    wc=p->ZP[KP1]-zp/p->DZP[KP];

    value =  lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::cint4a(double wx, double f0, double f1, double f2, double f3)
{
    double di0,di1,a0,a1,a2,a3,df;
    
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

