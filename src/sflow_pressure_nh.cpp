/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"sflow_pressure_nh.h"
#include"lexer.h"
#include"fdm2D.h"
#include"slice.h"
#include"sflow_curv.h"

#define HP (WL(i,j)>1.0e-20?WL(i,j):1.0e20)

int sflow_pressure_nh::is_active(lexer *p, fdm2D *b, int ii, int jj)
{
    i=ii;
    j=jj;
    return active(p,b);
}

int sflow_pressure_nh::active(lexer *p, fdm2D *b)
{
    // non-hydrostatic pressure is solved in wet, deep, non-breaking cells
    if(p->wet[IJ]==0 || p->deep[IJ]==0 || b->breaking(i,j)>0)
    return 0;
    
    if(p->wet[Im1J]==0 || p->wet[Ip1J]==0)
    return 0;
    
    if(p->j_dir==1 && (p->wet[IJm1]==0 || p->wet[IJp1]==0))
    return 0;
    
    return 1;
}

void sflow_pressure_nh::poisson(lexer*p, fdm2D* b, slice &WL, double cbed)
{
    n=0;

    // curvilinear grid: compact Laplacian of the orthogonal parts of the faces
    if(p->geo_curv==1)
    {
    double cE,cW,cN,cS;

	SLICELOOP4
	{
    sflow_curv::lapcoef(p,i,j,cE,cW,cN,cS);

	b->M.p[n] = WL(i,j)*(cE + cW + cN + cS)/p->W1 + 2.0*cbed/(HP*p->W1);
   	b->M.n[n] = -WL(i,j)*cE/p->W1;
	b->M.s[n] = -WL(i,j)*cW/p->W1;
	b->M.w[n] = -WL(i,j)*cN/p->W1;
	b->M.e[n] = -WL(i,j)*cS/p->W1;
	++n;
	}
    }

    if(p->geo_curv==0)
    SLICELOOP4
	{
	b->M.p[n] =  WL(i,j)/(p->W1*p->DXP[IP]*p->DXN[IP])
               + WL(i,j)/(p->W1*p->DXP[IM1]*p->DXN[IP])
               + WL(i,j)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir
               + WL(i,j)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir
               + 2.0*cbed/(HP*p->W1);
    
   	b->M.n[n] = -WL(i,j)/(p->W1*p->DXP[IP]*p->DXN[IP]);
	b->M.s[n] = -WL(i,j)/(p->W1*p->DXP[IM1]*p->DXN[IP]);
	b->M.w[n] = -WL(i,j)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir;
	b->M.e[n] = -WL(i,j)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir;
	++n;
	}
    
    n=0;
	SLICELOOP4
	{
        // inflow / outflow: q=0
		if(p->flagslice4[Im1J]<0 && p->IOSL[Im1J]==1)
		b->M.s[n] = 0.0;
        
        if(p->flagslice4[Ip1J]<0 && p->IOSL[Ip1J]==1)
		b->M.n[n] = 0.0;
        
        // walls: dq/dn=0
        if(p->flagslice4[Im1J]<0 && p->IOSL[Im1J]==0)
		{
		b->M.p[n] += b->M.s[n];
		b->M.s[n] = 0.0;
		}

        if(p->flagslice4[Ip1J]<0 && p->IOSL[Ip1J]==0)
		{
		b->M.p[n] += b->M.n[n];
		b->M.n[n] = 0.0;
		}

		if(p->flagslice4[IJm1]<0)
		{
		b->M.p[n] += b->M.e[n];
		b->M.e[n] = 0.0;
		}

		if(p->flagslice4[IJp1]<0)
		{
		b->M.p[n] += b->M.w[n];
		b->M.w[n] = 0.0;
		}
        
        // hydrostatic cells
        if(active(p,b)==0)
        {
        b->M.p[n] = 1.0;
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.w[n] = 0.0;
        b->M.e[n] = 0.0;
        b->rhsvec.V[n] = 0.0;
        b->press(i,j) = 0.0;
        }
	++n;
	}
}
