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

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void sediment_exner::timestep(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    double dx=1.0e8;
	maxvz=maxdh=0.0;
    
    
    SEDSLICELOOP
    {
        if(p->j_dir==1 && p->knoy>1)
        {
        dx = MIN(dx,p->DXN[IP]);
        dx = MIN(dx,p->DYN[JP]);
        }
        
        if(p->j_dir==0 || p->knoy==1)
        dx = MIN(dx,p->DXN[IP]);
    }
    
	SEDSLICELOOP
	maxvz = MAX(fabs(s->vz(i,j)),maxvz);

	maxvz=pgc->globalmax(maxvz);
	
    // 
    if(p->S29==0)
    {
	if(p->S15==0)
    p->dtsed=p->S17*MIN(p->S13, (p->S14*dx)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));

    if(p->S15==1)
    p->dtsed=MIN(p->S17*p->dt, (p->S14*dx)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S17*p->S13;
    
    if(p->S15==3)
    p->dtsed=p->S17*p->dt;
    }
    
    // ramp up
    if(p->S29==1)
    {
    p->dtsed = (p->S14*dx)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15);
    
    p->dtsed = MIN(p->dtsed,ramp_dt(p));
    }

    // bed celerity limit (S103>0)
    // Exner: dz/dt + S35/(1-n) * dqb/dx = 0  ->  celerity c = S35/(1-n) * dqb/dz
    // Roe-type estimate across each sediment face, only where |dz| > d50 so that
    // transport gradients over a flat bed do not produce spurious celerities.
    // dtsed <= S103 * dx/c
    if(p->S103>0.0)
    {
    double cdx=0.0;
    double dz,dq;
    const double dzmin = MAX(p->S20,1.0e-10);
    const double cfac = p->S35/(1.0-p->S24);

        SEDSLICELOOP
        {
            if(i+p->origin_i<p->gknox-1 && p->flagslice4[Ip1J]>0 && p->DFBED[Ip1J]>0)
            {
            dz = fabs(s->bedzh(i+1,j)-s->bedzh(i,j));
            dq = fabs(s->qb(i+1,j)-s->qb(i,j));

            if(dz>dzmin)
            cdx = MAX(cdx, cfac*dq/(dz*p->DXP[IP]));
            }

            if(p->j_dir==1 && p->gknoy>1)
            if(j+p->origin_j<p->gknoy-1 && p->flagslice4[IJp1]>0 && p->DFBED[IJp1]>0)
            {
            dz = fabs(s->bedzh(i,j+1)-s->bedzh(i,j));
            dq = fabs(s->qb(i,j+1)-s->qb(i,j));

            if(dz>dzmin)
            cdx = MAX(cdx, cfac*dq/(dz*p->DYP[JP]));
            }
        }

    cdx = pgc->globalmax(cdx);

    if(cdx>1.0e-20)
    p->dtsed = MIN(p->dtsed, p->S103/cdx);
    }

    p->dtsed=pgc->timesync(p->dtsed);
    
    //
    maxdh=p->dtsed*maxvz;
	
	if(p->mpirank==0)
	cout<<p->mpirank<<" max_vz: "<<setprecision(4)<<maxvz<<" max_dh: "<<setprecision(4)<<maxdh<<" dtsed: "<<setprecision(4)<<p->dtsed<<endl;
}

double sediment_exner::ramp_dt(lexer *p)
{
    double f=1.0;
    double dt=p->S13;

    if(p->sedtime>=p->S29_ts && p->sedtime<p->S29_te)
    {
    f = (p->sedtime-p->S29_ts)/(p->S29_te-p->S29_ts);
    
    dt = (1.0-f)*p->S29_dts + f*p->S29_dte;
    }

    if(p->sedtime<p->S29_ts)
    dt=p->S29_dts;
    
    if(p->sedtime>p->S29_te)
    dt=p->S29_dte;

    return dt;
}