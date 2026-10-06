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
#include"bedconc_VR.h"
#include"topo_relax.h"
#include"sediment_exnerdisc.h"
#include"sediment_fdm.h"
#include"sediment_mixture.h"

double sediment_exner::susp_ED(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    double val=0.0;
    
    if(p->S62==1 && p->count>p->S43)
    val = (-s->ws)*(s->cb(i,j) - s->cbe(i,j)); 
    
    // NHFLOW (S10 1): the suspended solver ran before bedconc and used the reference concentration
    // of the previous sediment step (cbn); erode the bed with the same value, otherwise every
    // wet/dry cycle of a swash cell leaves ws*dt*(cbe^n - cbe^(n-1)) unbalanced
    if(p->S62==1 && p->count>p->S43 && p->A10==5 && p->S10==1)
    val = (-s->ws)*(s->cb(i,j) - s->cbn(i,j)); 
    
    // NHFLOW: same exchange condition as the water column (nhflow_suspended_IM1::suspsource: k=0, wet, DF>0),
    // the bed cannot erode or deposit suspended sediment the water column never receives or loses
    if(p->A10==5)
    {
    k=0;
    if(p->wet[IJ]==0 || p->DF[IJK]<0 || p->DFBED[IJ]<0)
    val = 0.0;
    }
    
    // multi-fraction bed: exchange with the single suspended class weighted by F_k
    if(frac_k>=0)
    val *= (*s->pmix->F[frac_k])(i,j);
    
    // E: cbe
    // D: conc
    
    return val;
}

void sediment_exner::susp_qs(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    if(frac_k<0)
    SEDSLICELOOP
    s->qb(i,j) += s->qbs(i,j);
    
    if(frac_k>=0)
    SEDSLICELOOP
    s->qb(i,j) += (*s->pmix->F[frac_k])(i,j)*s->qbs(i,j);
    
    pgc->gcsl_start4(p,s->qb,1);
}