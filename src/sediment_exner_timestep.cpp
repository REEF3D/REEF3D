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
    
    // multi-fraction bed: fastest single fraction (set in start_mixture, 0 otherwise)
    maxvz = MAX(maxvz,maxvz_k);
    maxvz_k = 0.0;
	
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

    // flow time since the previous bed update (S 44 steps, S 46 seconds, ...)
    double dtflow = (simtime_last>=0.0 && p->simtime>simtime_last) ? p->simtime-simtime_last : p->dt;
    simtime_last = p->simtime;
    dtflow = pgc->timesync(dtflow);
    
    // S 15 4: morphological factor S 17 on the flow time since the last bed update
    // (S 15 3 uses the current flow step only: with S 44 > 1 the factor is S 17/S 44)
    if(p->S29==0 && p->S15==4)
    p->dtsed = p->S17*dtflow;
    
    // bed celerity (Roe-type estimate across each sediment face, only where |dz| > d50 so that
    // transport gradients over a flat bed do not produce spurious celerities)
    // Exner: dz/dt + S35/(1-n) * dqb/dx = 0  ->  celerity c = S35/(1-n) * dqb/dz
    // and the near-bed flow speed U = |(P,Q)|
    double cmax=0.0, cdx=0.0, umax=0.0;
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
            {
            cmax = MAX(cmax, cfac*dq/dz);
            cdx = MAX(cdx, cfac*dq/(dz*p->DXP[IP]));
            }
            }
            
            if(p->j_dir==1 && p->gknoy>1)
            if(j+p->origin_j<p->gknoy-1 && p->flagslice4[IJp1]>0 && p->DFBED[IJp1]>0)
            {
            dz = fabs(s->bedzh(i,j+1)-s->bedzh(i,j));
            dq = fabs(s->qb(i,j+1)-s->qb(i,j));
            if(dz>dzmin)
            {
            cmax = MAX(cmax, cfac*dq/dz);
            cdx = MAX(cdx, cfac*dq/(dz*p->DYP[JP]));
            }
            }
            
        umax = MAX(umax, sqrt(s->P(i,j)*s->P(i,j) + s->Q(i,j)*s->Q(i,j)));
        }
        
    cmax = pgc->globalmax(cmax);
    cdx = pgc->globalmax(cdx);
    umax = pgc->globalmax(umax);
    
    // bed celerity limit (S103>0): dtsed <= S103 * dx/c
    if(p->S103>0.0 && cdx>1.0e-20)
    p->dtsed = MIN(p->dtsed, p->S103/cdx);
    
    // morphological acceleration limit (S104>0): the bed forms move relative to the flow at most
    // at S104 times the near-bed flow speed, c*dtsed/dtflow <= S104*U; the error of the
    // morphological acceleration is of the order of this ratio (validation case 02)
    if(p->S104>0.0 && cmax>1.0e-20 && umax>1.0e-20)
    p->dtsed = MIN(p->dtsed, p->S104*umax*dtflow/cmax);
    
    p->dtsed=pgc->timesync(p->dtsed);
    
    // effective morphological factor and bed celerity / flow speed
    morfac = p->S35*p->dtsed/MAX(dtflow,1.0e-20);
    c_U = umax>1.0e-20 ? cmax*p->dtsed/(MAX(dtflow,1.0e-20)*umax) : 0.0;
    
    maxdh=p->dtsed*maxvz;
	
	if(p->mpirank==0)
	cout<<p->mpirank<<" max_vz: "<<setprecision(4)<<maxvz<<" max_dh: "<<setprecision(4)<<maxdh<<" dtsed: "<<setprecision(4)<<p->dtsed
        <<" morph. factor: "<<setprecision(4)<<morfac<<" c/U: "<<setprecision(3)<<c_U<<endl;
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