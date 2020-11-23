/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sediment_f::timestep(lexer *p, fdm *a, ghostcell *pgc)
{
    double vx,vy,vz;
    double xmax,ymax,zmax,vmax;

    xmax=ymax=zmax=vmax=0.0;
    
	SLICELOOP4
    {	
	topovel(p,a,pgc,vx,vy,vz);
	
    zmax=MAX(fabs(zmax),fabs(vz));
	}

    zmax=pgc->globalmax(zmax);
    
    if(p->S15==0)
    p->dtsed=MIN(p->S13, (p->S14*p->DXM)/(fabs(zmax)>1.0e-15?zmax:1.0e-15));

    if(p->S15==1)
    p->dtsed=MIN(p->dt, (p->S14*p->DXM)/(fabs(zmax)>1.0e-15?zmax:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S13;

    p->dtsed=pgc->timesync(p->dtsed);

    p->sediter=1;

    if(p->S17==1)
    p->sediter=MAX(int(p->S13/p->dtsed),1);

    p->sediter=MIN(p->sediter, p->S18);

    if(p->mpirank==0)
    cout<<endl<<"sediter: "<<p->sediter<<endl;

    cout<<p->mpirank<<" xmax: "<<setprecision(4)<<xmax<<" ymax: "<<setprecision(4)<<ymax<<" zmax: "<<setprecision(4)<<zmax<<" dtsed: "<<setprecision(4)<<p->dtsed<<endl;
}
