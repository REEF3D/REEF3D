/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"

void sflow_sediment_f::exner(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double uvel,vvel,u_abs;
	double signx,signy;
	double dqx,dqy;
    double maxvel;
	 
    SLICELOOP4
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{						
        pip=1;
        uvel=0.5*(b->P(i,j)+b->P(i-1,j));
        pip=0;

        pip=2;
        vvel=0.5*(b->Q(i,j)+b->Q(i,j-1));
        pip=0;
		
		u_abs = sqrt(uvel*uvel + vvel*vvel);
		signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
		signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;
	
		
	dqx=dqy=0.0;

    dqx = (b->qb(i+1,j)-b->qb(i-1,j))/(2.0*p->dx);
    dqy = (b->qb(i,j+1)-b->qb(i,j-1))/(2.0*p->dx);
		
	// Exner equation
    topovel(i,j) =  -(1.0/(1.0-p->S24))*(dqx*signx + dqy*signy); 
	}
    
    // timestep
    maxvel=0.0;
	SLICELOOP4
	maxvel = MAX(fabs(topovel(i,j)),maxvel);	
    
    if(p->S15==0)
    p->dtsed=MIN(p->S13, (p->S14*p->dx)/(fabs(maxvel)>1.0e-15?maxvel:1.0e-15));

    if(p->S15==1)
    p->dtsed=MIN(p->dt, (p->S14*p->dx)/(fabs(maxvel)>1.0e-15?maxvel:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S13;
    
    p->dtsed=pgc->timesync(p->dtsed);
    
     p->sedtime+=p->dtsed;
     
    // bedchange
    SLICELOOP4
    b->bed(i,j) += p->dtsed*topovel(i,j);
    
}