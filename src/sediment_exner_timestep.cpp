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
    double dx;
	maxvz=maxdh=0.0;
    
    
    SEDSLICELOOP
    {
    if(p->j_dir==1 && p->knoy>1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->j_dir==0 || p->knoy==1)
    dx = p->DXN[IP];
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
    }
    
    // ramp up
    if(p->S29==1)
    {
    p->dtsed = (p->S14*dx)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15);
    
    p->dtsed = MIN(p->dtsed,ramp_dt(p));
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