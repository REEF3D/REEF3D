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

#include"sflow_sediment_RK3.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"slice1.h"
#include"slice2.h"
#include"sflow_sediment_f.h"
 
sflow_sediment_RK3::sflow_sediment_RK3(lexer* p, fdm2D *b) : bedrk0(p),bedrk1(p),bedrk2(p)
{
    psed = new sflow_sediment_f(p,b);
}

sflow_sediment_RK3::~sflow_sediment_RK3()
{
}

void sflow_sediment_RK3::step1(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, double alpha)
{
    // bedrk0 ini
    SLICELOOP4
    bedrk0(i,j) = b->bed(i,j);
    
    // timestep calculation
    if(p->S15==0)
    p->dtsed=MIN(p->S13, (p->S14*p->dx)/(fabs(p->maxtopovel)>1.0e-15?p->maxtopovel:1.0e-15));

    if(p->S15==1)
    p->dtsed=MIN(p->dt, (p->S14*p->dx)/(fabs(p->maxtopovel)>1.0e-15?p->maxtopovel:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S13;
    
    p->dtsed=pgc->timesync(p->dtsed);
    
    p->sedtime+=p->dtsed;
    
    p->maxtopovel=0.0;
    
    // step 1
    psed->start(p,b,pgc,P,Q,b->topovel);
    
    // bedchange
    SLICELOOP4
    {
    bedrk1(i,j) = bedrk0(i,j) + p->dtsed*b->topovel(i,j);
    b->bed(i,j) = bedrk1(i,j);
    }
}

void sflow_sediment_RK3::step2(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, double alpha)
{
    // step 2
    psed->start(p,b,pgc,P,Q,b->topovel);
    
    // bedchange
    SLICELOOP4
    {
    bedrk2(i,j) = 0.75*bedrk0(i,j) + 0.25*bedrk1(i,j) + 0.25*p->dtsed*b->topovel(i,j); 
    b->bed(i,j) = bedrk2(i,j);
    }
}

void sflow_sediment_RK3::step3(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, double alpha)
{
    // step 3
    psed->start(p,b,pgc,P,Q,b->topovel);
    
    // bedchange
    SLICELOOP4
    b->bed(i,j) = (1.0/3.0)*bedrk0(i,j) + (2.0/3.0)*bedrk2(i,j) + (2.0/3.0)*p->dtsed*b->topovel(i,j); 
    
}

void sflow_sediment_RK3::step4(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, double alpha)
{
    
}