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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<math.h>

void particle_pls::correct(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow)
{
    corrected=0;
	
    inicorr(p,a,pgc);

    for(n=0;n<posactive;++n)
    if(posflag[n]>0)
	//if(pos[n][3]<=-pos[n][4])
    {
    ii=int(pos[n][0]/dx-0.5);
    jj=int(pos[n][1]/dx-0.5);
    kk=int(pos[n][2]/dx-0.5);

    parcorr(p,a,pflow,1.0,pos[n],ii,jj,kk,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii+1,jj,k,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii,jj+1,kk,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii+1,jj+1,kk,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii,jj,kk+1,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii+1,jj,kk+1,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii,jj+1,kk+1,ii,jj,kk);
    parcorr(p,a,pflow,1.0,pos[n],ii+1,jj+1,kk+1,ii,jj,kk);
    }

    for(n=0;n<negactive;++n)
    if(negflag[n]>0)
	//if(neg[n][3]>=neg[n][4])
    {
    ii=int(neg[n][0]/dx-0.5);
    jj=int(neg[n][1]/dx-0.5);
    kk=int(neg[n][2]/dx-0.5);

    parcorr(p,a,pflow,-1.0,neg[n],ii,jj,k,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii+1,jj,kk,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii,jj+1,kk,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii+1,jj+1,kk,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii,jj,kk+1,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii+1,jj,kk+1,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii,jj+1,kk+1,ii,jj,kk);
    parcorr(p,a,pflow,-1.0,neg[n],ii+1,jj+1,kk+1,ii,jj,kk);
    }

    finalcorr(p,a,pgc);
}

void particle_pls::parcorr(lexer *p,fdm* a,ioflow *pflow,double sign,double* f,int i1,int j1,int k1,int i0, int j0, int k0)
{
    int check=0;
	
    i=i1;
    j=j1;
    k=k1;
	
    check=boundcheck(p,a,i,j,k,1);
	
	if(check==1)
	check=pflow->iozonecheck(p,a);
	
    if(check==1)
    {
    xc=(double(i) + 0.5)*dx;
    yc=(double(j) + 0.5)*dx;
    zc=(double(k) + 0.5)*dx;

    xcell=(double(i0) + 0.5)*dx;
    ycell=(double(j0) + 0.5)*dx;
    zcell=(double(k0) + 0.5)*dx;

    normreg(a,i,j,k);
	
    length = sqrt(pow(f[0]-xc,2.0) + pow(f[1]-yc,2.0) + pow(f[2]-zc,2.0));

    gamma=1.0-(f[4]/(length>1.0e-15?length:1.0e20));

    scalar = (f[0]-xc)*nvec[0] + (f[1]-yc)*nvec[1] + (f[2]-zc)*nvec[2];

    xs = xc + gamma*nvec[0]*scalar;
    ys = yc + gamma*nvec[1]*scalar;
    zs = zc + gamma*nvec[2]*scalar;

		if(xs>=xcell-nu && xs<=xcell+dx+nu && ys>=ycell-nu && ys<=ycell+dx+nu && zs>=zcell-nu && zs<=zcell+dx+nu)
		{
			sp=-1.0;
			// +
			if(sign>0.0)
			{
				if(a->phi(i,j,k)<=0.0)
				phimax(i,j,k) = MAX(phimax(i,j,k), (f[4]-length));

				if(a->phi(i,j,k)>0.0)
				phimax(i,j,k) = MAX(phimax(i,j,k), (f[4]+length));
			}

			// -
			if(sign<0.0)
			{
				if(a->phi(i,j,k)>0.0)
				phimin(i,j,k) = MIN(phimin(i,j,k), sp*(f[4]-length));

				if(a->phi(i,j,k)<=0.0)
				phimin(i,j,k) = MIN(phimin(i,j,k), sp*(f[4]+length));
			}
		}
	}
}

void particle_pls::finalcorr(lexer* p, fdm* a, ghostcell* pgc)
{
		
	LOOP
    {	
		if(fabs(phimax(i,j,k))<fabs(phimin(i,j,k)))
		{
		a->phi(i,j,k)=phimax(i,j,k);
		
		corrected++;
		}

		// ---
		
		if(fabs(phimax(i,j,k))>fabs(phimin(i,j,k)))
		{
		a->phi(i,j,k)=phimin(i,j,k);
		
		corrected++;
		}
    }
}

void particle_pls::inicorr(lexer* p, fdm* a, ghostcell* pgc)
{
    LOOP
    {
    phimax(i,j,k)=a->phi(i,j,k);
    phimin(i,j,k)=a->phi(i,j,k);
	phiold(i,j,k)=a->phi(i,j,k);
    }
	
	pgc->start4(p,phimax,1);
	pgc->start4(p,phimin,1);
	pgc->start4(p,phiold,1);
}
