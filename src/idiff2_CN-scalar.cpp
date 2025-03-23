/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"idiff2_CN.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

void idiff2_CN::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field& b, field& visc, double sig, double alpha)
{

}

void idiff2_CN::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field& b, field& visc, field &eddyv, double sig, double alpha)
{
    starttime=pgc->timer();
    double visctermp, viscterms, visctermn, viscterme, visctermw, visctermb, visctermt;
    
    count=0;
	LOOP
	{
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=visc(i,j,k);

	
//   M
	
	visctermp  =    0.5*(visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + visc_ijk+ev_ijk/sig)/(p->DXN[IP]*p->DXP[IM1])
					+   0.5*(visc_ijk+ev_ijk/sig + visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)/(p->DXN[IP]*p->DXP[IP])
					+   0.5*(visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + visc_ijk+ev_ijk/sig)/(p->DYN[JP]*p->DYP[JM1])
					+   0.5*(visc_ijk+ev_ijk/sig + visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)/(p->DYN[JP]*p->DYP[JP])
					+   0.5*(visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + visc_ijk+ev_ijk/sig)/(p->DZN[KP]*p->DZP[KM1])
					+   0.5*(visc_ijk+ev_ijk/sig + visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)/(p->DZN[KP]*p->DZP[KP]);

         a->M.p[count]  = alphaCN*visctermp + 1.0/(alpha*p->dt);

         a->rhsvec.V[count] -= (1-alphaCN)*visctermp*b(i,j,k);
    
         a->rhsvec.V[count] += b(i,j,k)/(alpha*p->dt); 
	 
	 viscterms = -0.5*(visc_ijk+ev_ijk/sig + visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)/(p->DXN[IP]*p->DXP[IM1]);
	 visctermn = -0.5*(visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + visc_ijk+ev_ijk/sig)/(p->DXN[IP]*p->DXP[IP]);
	 viscterme = -0.5*(visc_ijk+ev_ijk/sig + visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)/(p->DYN[JP]*p->DYP[JM1]);
	 visctermw = -0.5*(visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + visc_ijk+ev_ijk/sig)/(p->DYN[JP]*p->DYP[JP]);
	 visctermb = -0.5*(visc_ijk+ev_ijk/sig + visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)/(p->DZN[KP]*p->DZP[KM1]);
	 visctermt = -0.5*(visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + visc_ijk+ev_ijk/sig)/(p->DZN[KP]*p->DZP[KP]);
	 
        a->M.s[count] = alphaCN*viscterms;
        a->M.n[count] = alphaCN*visctermn;
        a->M.e[count] = alphaCN*viscterme;
        a->M.w[count] = alphaCN*visctermw;
        a->M.b[count] = alphaCN*visctermb;
        a->M.t[count] = alphaCN*visctermt;

        a->rhsvec.V[count] -= (1-alphaCN)*(viscterms*b(i-1,j,k) + visctermn*b(i+1,j,k) + viscterme*b(i,j-1,k) + visctermw*b(i,j+1,k) + visctermb*b(i,j,k-1) + visctermt*b(i,j,k+1));

	 ++count;
	}
    
    n=0;
	LOOP
	{
		if(p->flag4[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*b(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag4[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*b(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*b(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*b(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*b(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*b(i,j,k+1);
		a->M.t[n] = 0.0;
		}

	++n;
	}
    
    psolv->start(p,a,pgc,b,a->rhsvec,4);
    time=pgc->timer()-starttime;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"scalar_diffiter: "<<p->solveriter<<"  scalar_difftime: "<<setprecision(3)<<time<<endl;
}



void idiff2_CN::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field& diff, field& b, field& visc, field &eddyv, double sig, double alpha)
{

    starttime=pgc->timer();
    double visctermp, viscterms, visctermn, viscterme, visctermw, visctermb, visctermt;

    count=0;
	LOOP
	{
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=visc(i,j,k);

	
//   M
	visctermp  =    0.5*(visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + visc_ijk+ev_ijk/sig)/(p->DXN[IP]*p->DXP[IM1])
                                        +   0.5*(visc_ijk+ev_ijk/sig + visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)/(p->DXN[IP]*p->DXP[IP])
                                        +   0.5*(visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + visc_ijk+ev_ijk/sig)/(p->DYN[JP]*p->DYP[JM1])
                                        +   0.5*(visc_ijk+ev_ijk/sig + visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)/(p->DYN[JP]*p->DYP[JP])
                                        +   0.5*(visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + visc_ijk+ev_ijk/sig)/(p->DZN[KP]*p->DZP[KM1])
                                        +   0.5*(visc_ijk+ev_ijk/sig + visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)/(p->DZN[KP]*p->DZP[KP]);

         a->M.p[count]  = alphaCN*visctermp + 1.0/(alpha*p->dt);

         a->rhsvec.V[count] -= (1-alphaCN)*visctermp*b(i,j,k);

         a->rhsvec.V[count] += b(i,j,k)/(alpha*p->dt);

         viscterms = -0.5*(visc_ijk+ev_ijk/sig + visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)/(p->DXN[IP]*p->DXP[IM1]);
         visctermn = -0.5*(visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + visc_ijk+ev_ijk/sig)/(p->DXN[IP]*p->DXP[IP]);
         viscterme = -0.5*(visc_ijk+ev_ijk/sig + visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)/(p->DYN[JP]*p->DYP[JM1]);
         visctermw = -0.5*(visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + visc_ijk+ev_ijk/sig)/(p->DYN[JP]*p->DYP[JP]);
         visctermb = -0.5*(visc_ijk+ev_ijk/sig + visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)/(p->DZN[KP]*p->DZP[KM1]);
         visctermt = -0.5*(visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + visc_ijk+ev_ijk/sig)/(p->DZN[KP]*p->DZP[KP]);

        a->M.s[count] = alphaCN*viscterms;
        a->M.n[count] = alphaCN*visctermn;
        a->M.e[count] = alphaCN*viscterme;
        a->M.w[count] = alphaCN*visctermw;
        a->M.b[count] = alphaCN*visctermb;
        a->M.t[count] = alphaCN*visctermt;

        a->rhsvec.V[count] -= (1-alphaCN)*(viscterms*b(i-1,j,k) + visctermn*b(i+1,j,k) + viscterme*b(i,j-1,k) + visctermw*b(i,j+1,k) + visctermb*b(i,j,k-1) + visctermt*b(i,j,k+1));
	 
	++count;
	}
    
    n=0;
	LOOP
	{
		if(p->flag4[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*b(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag4[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*b(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*b(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*b(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*b(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*b(i,j,k+1);
		a->M.t[n] = 0.0;
		}

	++n;
	}
    
    psolv->start(p,a,pgc,diff,a->rhsvec,4);
    time=pgc->timer()-starttime;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"scalar_diffiter: "<<p->solveriter<<"  scalar_difftime: "<<setprecision(3)<<time<<endl;
}
