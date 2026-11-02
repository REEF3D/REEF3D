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
#include"sediment_curv.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void sediment_exner::timestep(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    // curvilinear grid (G 50): sediment_curv.h
    if(p->geo_curv==1)
    {
    timestep_curv(p,pgc,s);
    return;
    }

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
    
    // mesh refinement, multi-block grids: the other grids (sediment_exner_amr.cpp)
    if(amr_maxvz>0.0)
    maxvz = MAX(maxvz,amr_maxvz);
    if(amr_mindx>0.0)
    dx = MIN(dx,amr_mindx);
    
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
    const double cfac = p->S35/(1.0-p->S24);
    celerity(p,s,nullptr,cmax,cdx,umax);

    // mesh refinement, multi-block grids: the other grids (amr_celerity, sflow_amr_sed.cpp, nhflow_amr_sed.cpp), so
    // that S 103 and S 104 give one dtsed for the whole hierarchy (before P9 only the grid of level 0 counted)
    cmax = MAX(cmax,amr_cmax);
    cdx = MAX(cdx,amr_cdx);
    umax = MAX(umax,amr_umax);

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

    // S 86: the slope flux -qb/f dz/dn is an explicit diffusion of the bed: dtsed <= 0.4/(S35/(1-n) max(qb/f (1/dx^2+1/dy^2)))
    if(p->S86>0)
    {
    const double rate = pgc->globalmax(MAX(bend_rate,amr_bendrate));
    if(rate>1.0e-20)
    p->dtsed = MIN(p->dtsed, 0.4/(cfac*rate));
    }
    
    p->dtsed=pgc->timesync(p->dtsed);
    
    // effective morphological factor and bed celerity / flow speed
    morfac = p->S35*p->dtsed/MAX(dtflow,1.0e-20);
    c_U = umax>1.0e-20 ? cmax*p->dtsed/(MAX(dtflow,1.0e-20)*umax) : 0.0;
    cu_warn(p);
    
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

// ---- curvilinear grid (G 50): the functions above with the geometry of sediment_curv.h

void sediment_exner::timestep_curv(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    double dx=1.0e8;
	maxvz=maxdh=0.0;
    
    
    SEDSLICELOOP
    {
        if(p->geo_curv==1)
        {
        dx = MIN(dx,sediment_curv::width(p,1,i,j));
        dx = MIN(dx,sediment_curv::width(p,2,i,j));
        }
        else
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
    
    // mesh refinement, multi-block grids: the other grids (sediment_exner_amr.cpp)
    if(amr_maxvz>0.0)
    maxvz = MAX(maxvz,amr_maxvz);
    if(amr_mindx>0.0)
    dx = MIN(dx,amr_mindx);
    
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
    const double cfac = p->S35/(1.0-p->S24);
    celerity(p,s,nullptr,cmax,cdx,umax);

    // mesh refinement, multi-block grids: the other grids (amr_celerity, sflow_amr_sed.cpp, nhflow_amr_sed.cpp), so
    // that S 103 and S 104 give one dtsed for the whole hierarchy (before P9 only the grid of level 0 counted)
    cmax = MAX(cmax,amr_cmax);
    cdx = MAX(cdx,amr_cdx);
    umax = MAX(umax,amr_umax);

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

    // S 86: the slope flux -qb/f dz/dn is an explicit diffusion of the bed: dtsed <= 0.4/(S35/(1-n) max(qb/f (1/dx^2+1/dy^2)))
    if(p->S86>0)
    {
    const double rate = pgc->globalmax(MAX(bend_rate,amr_bendrate));
    if(rate>1.0e-20)
    p->dtsed = MIN(p->dtsed, 0.4/(cfac*rate));
    }
    
    p->dtsed=pgc->timesync(p->dtsed);
    
    // effective morphological factor and bed celerity / flow speed
    morfac = p->S35*p->dtsed/MAX(dtflow,1.0e-20);
    c_U = umax>1.0e-20 ? cmax*p->dtsed/(MAX(dtflow,1.0e-20)*umax) : 0.0;
    cu_warn(p);
    
    maxdh=p->dtsed*maxvz;
	
	if(p->mpirank==0)
	cout<<p->mpirank<<" max_vz: "<<setprecision(4)<<maxvz<<" max_dh: "<<setprecision(4)<<maxdh<<" dtsed: "<<setprecision(4)<<p->dtsed
        <<" morph. factor: "<<setprecision(4)<<morfac<<" c/U: "<<setprecision(3)<<c_U<<endl;
}

// the bed celerity of the grid (rank-local maxima; Roe-type estimate across each face of two sediment cells, only where
// |dz| > d50 so that transport gradients over a flat bed give no celerity): c = S35/(1-n) |dqb|/|dz|, c/dx, and the
// near-bed flow speed |(P,Q)|.  Every face of an own cell to an erodible cell of the domain (bend_open: the cells around a
// patch or block from the flags of their source, outside the domain closed), from both sides, so that the faces of the
// interfaces of a multi-block grid count as in one block; mask (mesh refinement): the cells with mask > 0.5 only
void sediment_exner::celerity(lexer *p, sediment_fdm *s, slice *mask, double &cmax, double &cdx, double &umax)
{
    const double dzmin = MAX(p->S20,1.0e-10);
    const double cfac = p->S35/(1.0-p->S24);
    const bool ydir = (p->j_dir==1 && p->gknoy>1);

    auto face = [&](int a, int b, double dl)
    {
        if(!bend_open(p,a,b) || (mask!=nullptr && (*mask)(a,b)<0.5))
        return;
        const double dz = fabs(s->bedzh(a,b)-s->bedzh(i,j));
        const double dq = fabs(s->qb(a,b)-s->qb(i,j));
        if(dz>dzmin)
        {
        cmax = MAX(cmax, cfac*dq/dz);
        cdx = MAX(cdx, cfac*dq/(dz*dl));
        }
    };

    SEDSLICELOOP
    if(mask==nullptr || (*mask)(i,j)>0.5)
    {
        face(i+1,j,p->geo_curv==1 ? sediment_curv::dist(p,i,j,1,0) : p->DXP[IP]);
        face(i-1,j,p->geo_curv==1 ? sediment_curv::dist(p,i-1,j,1,0) : p->DXP[IM1]);

        if(ydir)
        {
        face(i,j+1,p->geo_curv==1 ? sediment_curv::dist(p,i,j,0,1) : p->DYP[JP]);
        face(i,j-1,p->geo_curv==1 ? sediment_curv::dist(p,i,j-1,0,1) : p->DYP[JM1]);
        }

        umax = MAX(umax, sqrt(s->P(i,j)*s->P(i,j) + s->Q(i,j)*s->Q(i,j)));
    }
}

// the morphological acceleration beyond its range: the bed moves faster than the flow (c/U > 2, c the bed celerity times
// the morphological factor), once per run and only without the limits S 103 / S 104.  In the LFM bend (t9) MF 100 gives
// c/U 3 - 38 and bed waves that run upstream through the inflow reach in every variant (S 95 1, S 95 2, S 64); MF 10
// c/U 0.15.  The log only, the results stay the same.
void sediment_exner::cu_warn(lexer *p)
{
    if(cu_warned || c_U<=2.0 || p->S103>0.0 || p->S104>0.0)
    return;

    cu_warned = 1;
    if(p->mpirank==0)
    cout<<"!!! S 35: the bed moves faster than the flow (bed celerity x morphological factor / flow speed c/U = "
        <<setprecision(3)<<c_U<<" > 2): the morphological acceleration is beyond its range; lower S 35 or limit dtsed with "
        <<"S 104 (e.g. 0.5) or S 103 !!!"<<endl;
}
