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

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"sediment_curv.h"

void sediment_f::prep_cfd(lexer *p, fdm *a,ghostcell *pgc)
{    
    // vel prep --------
    pgc->start1(p,a->u,14);
	pgc->start2(p,a->v,15);
	pgc->start3(p,a->w,16);
    
    // find bedk -------
    fill_bedk(p,a,pgc);
    
    fill_PQ_cfd(p,a,pgc);
    
    waterlevel(p,a,pgc);
}

void sediment_f::prep_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc)
{    
    fill_PQ_nhflow(p,d,pgc);
    
    // S 37 2: bed shear stress (S 16 1) and transport direction (P, Q) from the phase-lead velocity
    if(p->S37==2 && p->geo_curv==0)
    phase_lead_nhflow(p,d,pgc);
}

/*--------------------------------------------------------------------
Thin water (NHFLOW: swash, wetting and drying fronts, shoals). The bed shear stress of a layer of a few grain
roughnesses is large (c_f of order 0.1) with any friction law, and the velocity of a front is the least reliable
one of the flow, so that a transport formula returns its largest rates where the least water flows.

S 131: the flow cannot carry more sediment than a share of its own volume flux,
           qb <= S131 h |u|,
       h the water depth and u the velocity of the bottom layer. The limit goes to zero with the depth and with the
       velocity and is far above the transport of deeper water (h = 1 m, 1 m/s, S131 0.05: 0.05 m^2/s).
S 132: no bed load and no pick-up of suspended load below the depth S132, linear to the full values at 2 S132;
       deposition goes on. The fluxes are scaled, not the bed change: the bed volume is conserved.
The bed load of a graded bed (S 51) is not limited (its fractions are set in sediment_mixture), the pick-up is.
conc 1: the reference concentration cbe as well (once per sediment step).
--------------------------------------------------------------------*/
void sediment_f::thin_water(lexer *p, ghostcell *pgc, int conc)
{
    if(p->S131<=0.0 && p->S132<=0.0)
    return;
    
    double h,u,v,f;
    
    SEDSLICELOOP
    {
    h = MAX(s->waterlevel(i,j),0.0);
    
    f = 1.0;
    
    if(p->S132>0.0)
    f = MAX(0.0, MIN(1.0, h/p->S132 - 1.0));
    
        if(s->pmix==nullptr)
        {
        s->qbe(i,j) *= f;
        
            if(p->S131>0.0)
            {
                if(p->geo_curv==1)
                sediment_curv::cellvel(p,s,i,j,u,v);
                else
                {
                u = 0.5*(s->P(i,j)+s->P(i-1,j));
                v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
                }
            
            s->qbe(i,j) = MIN(s->qbe(i,j), p->S131*h*sqrt(u*u + v*v));
            }
        }
    
    if(conc==1)
    s->cbe(i,j) *= f;
    }
    
    pgc->gcsl_start4(p,s->qbe,1);
}

void sediment_f::prep_sflow(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    fill_PQ_sflow(p,b,pgc,P,Q);
}