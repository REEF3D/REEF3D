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

void sediment_f::active_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    
    SLICEBASELOOP
    s->active(i,j)=0;
    
    ILOOP
    JLOOP
    {
        KWLOOP
        PSOLIDCHECK
        {
            
        if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
        s->active(i,j)=1;
        }
    }
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1); 
    
    // assign gcsldfbed entries
    SLICEBASELOOP
    {
    k = s->bedk(i,j);
    
    if(p->DF[IJK]>0)
    p->DFBED[IJ]=1;
    
    if(p->DF[IJK]<0)
    p->DFBED[IJ]=-1;
    }
    
    ALOOP
    a->test(i,j,k) = p->DFBED[IJ];
}

void sediment_f::active_ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    SLICEBASELOOP
    s->active(i,j)=1;
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1);
    
    
    // assign gcsldfbed entries
    SLICEBASELOOP
    {
    k = s->bedk(i,j);
    
    if(p->DF[IJK]>0)
    p->DFBED[IJ]=1;
    
    if(p->DF[IJK]<0)
    p->DFBED[IJ]=-1;
    }
}

void sediment_f::active_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    SLICEBASELOOP
    s->active(i,j)=1;
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1);
    
    
    // assign gcsldfbed entries
    k=0;
    SLICEBASELOOP
    {
    if(p->DF[IJK]>0)
    p->DFBED[IJ]=1;
    
    if(p->DF[IJK]<0)
    p->DFBED[IJ]=-1;
    }
}

void sediment_f::active_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    SLICEBASELOOP
    s->active(i,j)=1;
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1);
    
    
    // assign gcsldfbed entries
    k=0;
    SLICEBASELOOP
    {
    if(p->DF[IJK]>0)
    p->DFBED[IJ]=1;
    
    if(p->DF[IJK]<0)
    p->DFBED[IJ]=-1;
    }
}

void sediment_f::active_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICEBASELOOP
    s->active(i,j)=1;
    
    SLICEBASELOOP
    if(b->solidbed(i,j) >= s->bedzh(i,j))
    {
    s->active(i,j)=0;
    }
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1);
    
    
    // assign gcsldfeta entries
    k=p->knoz-1;
    SLICEBASELOOP
    {
    if(p->DF[IJK]>0)
    p->DFBED[IJ]=1;
    
    if(p->DF[IJK]<0)
    p->DFBED[IJ]=-1;
    }
}

void sediment_f::active_ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICEBASELOOP
    s->active(i,j)=1;
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1);
    
    
    // assign gcsldfeta entries
    k=p->knoz-1;
    SLICEBASELOOP
    {
    if(p->DF[IJK]>0)
    p->DFBED[IJ]=1;
    
    if(p->DF[IJK]<0)
    p->DFBED[IJ]=-1;
    }
}

void sediment_f::active_zone(lexer *p, ghostcell *pgc)
{
    SLICEBASELOOP
    for(int n=0;n<p->S74;++n)
    if(p->XP[IP]>p->S74_xs[n] && p->XP[IP]<p->S74_xe[n] && p->YP[JP]>p->S74_ys[n] && p->YP[JP]<p->S74_ye[n])
    s->active(i,j)=0;

}