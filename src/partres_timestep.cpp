/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void partres::timestep(lexer *p, ghostcell *pgc)
{
    double maxVelU=.00,maxVelV=0.0,maxVelW=0.0;
    double maxvz=0.0;
    
    for(size_t n=0;n<P.index;n++)
    {
        if(P.Flag[n]>=0)
        {
                maxVelU = MAX(maxVelU,fabs(P.U[n]));
                maxVelV = MAX(maxVelV,fabs(P.V[n]));
                maxVelW = MAX(maxVelW,fabs(P.W[n]));
        }
    }
    
    
    maxvz = MAX(maxVelU,maxVelV);
    maxvz = MAX(maxvz,maxVelV);
    
    maxvz = pgc->globalmax(maxvz);
    
    maxVelU = pgc->globalmax(maxVelU);
    maxVelV = pgc->globalmax(maxVelV);
    maxVelW = pgc->globalmax(maxVelW);
    
    if(timestep_ini==0)
    {
    maxvz = 1000.0;
    timestep_ini=1;
    }
    
    
    if(p->S15==0)
    p->dtsed=MIN(p->S13, (p->S14*p->DXM)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));

    if(p->S15==1)
    p->dtsed=MIN(p->dt, (p->S14*p->DXM)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S13;

    p->dtsed=pgc->timesync(p->dtsed);
    
    //
	
	if(p->mpirank==0)
    {
	cout<<p->mpirank<<" maxvz: "<<setprecision(4)<<maxvz<<" dtsed: "<<setprecision(4)<<p->dtsed<<endl;
    cout<<"Up_max: "<<maxVelU<<endl;
    cout<<"Vp_max: "<<maxVelV<<endl;
    cout<<"Wp_max: "<<maxVelW<<endl;
    }
    
    
}