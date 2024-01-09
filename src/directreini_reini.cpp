/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"directreini.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void directreini::reini(lexer *p,fdm* a, ghostcell *pgc, field& b, fieldint& nodeflag, fieldint& vertice)
{
    double xp[2],yp[2],zp[2];
    double nx,ny,nz;
    double nl;
    double d,t0,dist;
    double x0,y0,z0;
    double det0,det1,det2,det3,det4;
    int check,qn;

    for(n=0;n<numtri; ++n)
    if(numfac[n]>0)
	for(q=0;q<4;++q)	
    {
    ls[tri[n][q]] = sg(ls[tri[n][q]])*1.0e20;
	
	lsvert[tri[n][q]] = sg(lsvert[tri[n][q]])*1.0e20;
	
	lsfac[tri[n][q]] = sg(lsfac[tri[n][q]])*1.0e20;

	reiniflag[tri[n][q]]=0;	
    }

	int checker=0;
    for(n=0;n<numtri; ++n)
    if(numfac[n]>0)
    for(q=0;q<4;++q)
    {

    xp[0]=ccpt[facet[confac[n]][1]][0]-ccpt[facet[confac[n]][0]][0];
    yp[0]=ccpt[facet[confac[n]][1]][1]-ccpt[facet[confac[n]][0]][1];
    zp[0]=ccpt[facet[confac[n]][1]][2]-ccpt[facet[confac[n]][0]][2];

    xp[1]=ccpt[facet[confac[n]][2]][0]-ccpt[facet[confac[n]][0]][0];
    yp[1]=ccpt[facet[confac[n]][2]][1]-ccpt[facet[confac[n]][0]][1];
    zp[1]=ccpt[facet[confac[n]][2]][2]-ccpt[facet[confac[n]][0]][2];

    nx=yp[0]*zp[1] - zp[0]*yp[1];
    ny=zp[0]*xp[1] - xp[0]*zp[1];
    nz=xp[0]*yp[1] - yp[0]*xp[1];
	
	nl = sqrt(nx*nx + ny*ny + nz*nz);


    d = nx*ccpt[facet[confac[n]][0]][0] + ny*ccpt[facet[confac[n]][0]][1] + nz*ccpt[facet[confac[n]][0]][2];

    t0 = (d - (nx*pt[tri[n][q]][0] + ny*pt[tri[n][q]][1] + nz*pt[tri[n][q]][2]))
            /(nx*nx + ny*ny + nz*nz);

    x0 = pt[tri[n][q]][0] + t0*nx;
    y0 = pt[tri[n][q]][1] + t0*ny;
    z0 = pt[tri[n][q]][2] + t0*nz;


    det0 = determinant(pt[tri[n][0]][0], pt[tri[n][0]][1], pt[tri[n][0]][2], 1.0,
                       pt[tri[n][1]][0], pt[tri[n][1]][1], pt[tri[n][1]][2], 1.0,
                       pt[tri[n][2]][0], pt[tri[n][2]][1], pt[tri[n][2]][2], 1.0,
                       pt[tri[n][3]][0], pt[tri[n][3]][1], pt[tri[n][3]][2], 1.0);

    det1 = determinant(x0, y0, z0, 1.0,
                       pt[tri[n][1]][0], pt[tri[n][1]][1], pt[tri[n][1]][2], 1.0,
                       pt[tri[n][2]][0], pt[tri[n][2]][1], pt[tri[n][2]][2], 1.0,
                       pt[tri[n][3]][0], pt[tri[n][3]][1], pt[tri[n][3]][2], 1.0);

    det2 = determinant(pt[tri[n][0]][0], pt[tri[n][0]][1], pt[tri[n][0]][2], 1.0,
                       x0, y0, z0, 1.0,
                       pt[tri[n][2]][0], pt[tri[n][2]][1], pt[tri[n][2]][2], 1.0,
                       pt[tri[n][3]][0], pt[tri[n][3]][1], pt[tri[n][3]][2], 1.0);

    det3 = determinant(pt[tri[n][0]][0], pt[tri[n][0]][1], pt[tri[n][0]][2], 1.0,
                       pt[tri[n][1]][0], pt[tri[n][1]][1], pt[tri[n][1]][2], 1.0,
                       x0, y0, z0, 1.0,
                       pt[tri[n][3]][0], pt[tri[n][3]][1], pt[tri[n][3]][2], 1.0);

    det4 = determinant(pt[tri[n][0]][0], pt[tri[n][0]][1], pt[tri[n][0]][2], 1.0,
                       pt[tri[n][1]][0], pt[tri[n][1]][1], pt[tri[n][1]][2], 1.0,
                       pt[tri[n][2]][0], pt[tri[n][2]][1], pt[tri[n][2]][2], 1.0,
                       x0, y0, z0, 1.0);


    check=0;

    if(det0>=0.0 && det1>=0.0 && det2>=0.0 && det3>=0.0 && det4>=0.0)
    check=1;

    if(det0<=0.0 && det1<=0.0 && det2<=0.0 && det3<=0.0 && det4<=0.0)
    check=1;
	
	if(check==1)
	++checker;

        if(check==1)
		{
        dist = sqrt(pow(x0-pt[tri[n][q]][0], 2.0) + pow(y0-pt[tri[n][q]][1], 2.0) + pow(z0-pt[tri[n][q]][2], 2.0));
		
		++reiniflag[tri[n][q]];
		lsfac[tri[n][q]] = sg(lsfac[tri[n][q]])*MIN(fabs(lsfac[tri[n][q]]), fabs(dist));
		}
		

        if(check==0)
		if(reiniflag[tri[n][q]]==0)
        {
        dist=1.0e20;
        for(qn=0;qn<numfac[n];++qn)
        dist = MIN(dist,
               sqrt(pow(ccpt[facet[confac[n]][qn]][0]-pt[tri[n][q]][0], 2.0)
                  + pow(ccpt[facet[confac[n]][qn]][1]-pt[tri[n][q]][1], 2.0)
                  + pow(ccpt[facet[confac[n]][qn]][2]-pt[tri[n][q]][2], 2.0)));
				  
		lsvert[tri[n][q]] = sg(lsvert[tri[n][q]])*MIN(fabs(lsvert[tri[n][q]]), fabs(dist));
        }
	
	//ls[tri[n][q]] = sg(ls[tri[n][q]])*MIN(fabs(ls[tri[n][q]]), fabs(dist));
    
    }
	
	for(n=0;n<numtri;++n)
    for(q=0;q<4;++q)
    if(numfac[n]>0)
	{
		if(reiniflag[tri[n][q]]>=1)
		{
		ls[tri[n][q]]=lsfac[tri[n][q]];	
		}
		
		if(reiniflag[tri[n][q]]==0)
		ls[tri[n][q]]=fabs(lsvert[tri[n][q]])<fabs(lsfac[tri[n][q]])?lsvert[tri[n][q]]:lsfac[tri[n][q]];
		
	}
	
	checker=pgc->globalisum(checker);

	if(p->mpirank==0)
	cout<<"DirectReini Corrections: "<<checker<<endl;
	/*
    for(n=0;n<numtri;++n)
    for(q=0;q<4;++q)
    if(numfac[n]>0)
    {
    i=ijk[tri[n][q]][0];
    j=ijk[tri[n][q]][1];
    k=ijk[tri[n][q]][2];
    b(i,j,k)=ls[tri[n][q]];
    }
	*/
	for(n=0;n<numvert;++n)
    {
    i=ijk[n][0];
    j=ijk[n][1];
    k=ijk[n][2];
    b(i,j,k)=ls[n];
    }
}

double directreini::sg(double d)
{
    double s=0.0;

    if(d<0.0)
    s=-1.0;

    if(d>=0.0)
    s=1.0;

    return s;
}

double directreini::determinant(double m00,double m01,double m02,double m03,
                                double m10,double m11,double m12,double m13,
                                double m20,double m21,double m22,double m23,
                                double m30,double m31,double m32,double m33)
{
      double val=0.0;

      val =
      m03*m12*m21*m30 - m02*m13*m21*m30 - m03*m11*m22*m30 + m01*m13*m22*m30 +
      m02*m11*m23*m30 - m01*m12*m23*m30 - m03*m12*m20*m31 + m02*m13*m20*m31 +
      m03*m10*m22*m31 - m00*m13*m22*m31 - m02*m10*m23*m31 + m00*m12*m23*m31 +
      m03*m11*m20*m32 - m01*m13*m20*m32 - m03*m10*m21*m32 + m00*m13*m21*m32 +
      m01*m10*m23*m32 - m00*m11*m23*m32 - m02*m11*m20*m33 + m01*m12*m20*m33 +
      m02*m10*m21*m33 - m00*m12*m21*m33 - m01*m10*m22*m33 + m00*m11*m22*m33;

      return val;
}
