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
#include<fstream>

void directreini::debug(lexer* p,fdm* a)
{

	double xp[6],yp[6],zp[6];
    double nx,ny,nz;
    double nl;
    double d,t0,d_old;
    double x0,y0,z0;

    char name[100];
    sprintf(name,"directreini_debug-%i.dat",p->mpirank+1);
	ofstream result;
	result.open(name);
    count=0;
/*
    for(n=0;n<numtri; ++n)
    {
    result<<"n: "<<n<<" tri: "<<tri[n][0]<<" "<<tri[n][1]<<" "<<tri[n][2]<<" "<<tri[n][3]<<" confac: "<<confac[n]<<" numfac: "<<numfac[n];
    result<<"\t i: "<<ijk[tri[n][0]][0]<<" j: "<<ijk[tri[n][0]][1]<<" k: "<<ijk[tri[n][0]][2];
    result<<"\t x: "<<pt[tri[n][0]][0]<<" y: "<<pt[tri[n][0]][1]<<" z: "<<pt[tri[n][0]][2];
    result<<"\t ls1: "<<ls[tri[n][0]]<<" ls2: "<<ls[tri[n][1]]<<" ls3: "<<ls[tri[n][2]]<<" ls4: "<<ls[tri[n][3]]<<endl;
    }
*/
    for(n=0;n<numtri; ++n)
    {
    result<<"\t n: "<<n<<endl;
    for(q=0;q<4;++q)
    {

    result<<"\t x: "<<pt[tri[n][q]][0]<<" y: "<<pt[tri[n][q]][1]<<" z: "<<pt[tri[n][q]][2]<<" ls: "<<ls[tri[n][q]]<<endl;
    }
    for(int q=0;q<numfac[n];++q)
    {
    result<<ccpt[facet[confac[n]][q]][0]<<"  "<<ccpt[facet[confac[n]][q]][1]<<"  "<<ccpt[facet[confac[n]][q]][2]<<"  "<<endl;

    xp[0]=ccpt[facet[confac[n]][1]][0]-ccpt[facet[confac[n]][0]][0];
    yp[0]=ccpt[facet[confac[n]][1]][1]-ccpt[facet[confac[n]][0]][1];
    zp[0]=ccpt[facet[confac[n]][1]][2]-ccpt[facet[confac[n]][0]][2];

    xp[1]=ccpt[facet[confac[n]][2]][0]-ccpt[facet[confac[n]][0]][0];
    yp[1]=ccpt[facet[confac[n]][2]][1]-ccpt[facet[confac[n]][0]][1];
    zp[1]=ccpt[facet[confac[n]][2]][2]-ccpt[facet[confac[n]][0]][2];

    nx=yp[0]*zp[1] - zp[0]*yp[1];
    ny=zp[0]*xp[1] - xp[0]*zp[1];
    nz=xp[0]*yp[1] - yp[0]*xp[1];

    d = nx*ccpt[facet[confac[n]][0]][0] + ny*ccpt[facet[confac[n]][0]][1] + nz*ccpt[facet[confac[n]][0]][2];

    nl = sqrt(nx*nx + ny*ny + nz*nz);

    nl=fabs(nl)>1.0e-10?nl:1.0e10;
    d_old=d;

    nx/=nl;
    ny/=nl;
    nz/=nl;
    d/=nl;

    if(d<0.0)
    {
    nx*=-1.0;
    ny*=-1.0;
    nz*=-1.0;
    d*=-1.0;
    }

    t0 = (d - (nx*pt[tri[n][q]][0] + ny*pt[tri[n][q]][1] + nz*pt[tri[n][q]][2]))
            /(nx*nx + ny*ny + nz*nz);

    x0 = pt[tri[n][q]][0] + t0*nx;
    y0 = pt[tri[n][q]][1] + t0*ny;
    z0 = pt[tri[n][q]][2] + t0*nz;
    
    result<<"x: "<<x0<<" "<<y0<<" "<<z0<<endl;
    result<<"nx: "<<nx<<" ny: "<<ny<<" nz: "<<nz<<" d: "<<d<<" nl: "<<nl<<endl;
    result<<"b: "<<xp[0]<<" "<<yp[0]<<" "<<zp[0]<<"   c: "<<xp[1]<<" "<<yp[1]<<" "<<zp[1]<<endl;
    result<<"d_old: "<<d_old<<" d: "<<d<<" . "<<nl<<endl<<endl;
    }
    result<<" -------"<<endl;
    }

}
