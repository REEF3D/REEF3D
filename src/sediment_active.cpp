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
    
    LOOP
    a->test(i,j,k) = p->DFBED[IJ];
    
    dfbed_comms(p,pgc);
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
    
    dfbed_comms(p,pgc);
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
    
    dfbed_comms(p,pgc);
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
    
    dfbed_comms(p,pgc);
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
    
    
    // sediment cells: SFLOW has no 3D grid flags (DF), use the 2D cell flag
    SLICEBASELOOP
    {
    if(p->flagslice4[IJ]>0)
    p->DFBED[IJ]=1;
    
    if(p->flagslice4[IJ]<0)
    p->DFBED[IJ]=-1;
    }
    
    for(int q : dfbed_open)
    p->DFBED[q]=1;
    
    dfbed_inflow(p);
    
    dfbed_comms(p,pgc);
}

void sediment_f::active_ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // flagini() (3D grids) is not called for SFLOW: allocate the sediment cell flag here
    if(p->DFBED==nullptr)
    {
    p->Iarray(p->DFBED,p->imax*p->jmax);
    
    for(i=-p->margin; i<p->knox+p->margin; ++i)
    for(j=-p->margin; j<p->knoy+p->margin; ++j)
    p->DFBED[(i-p->imin)*p->jmax + j-p->jmin] = 1;
    }
    
    SLICEBASELOOP
    s->active(i,j)=1;
    
    active_zone(p,pgc);
    
    pgc->gcsl_start4int(p,s->active,1);
    
    
    // sediment cells: SFLOW has no 3D grid flags (DF), use the 2D cell flag
    SLICEBASELOOP
    {
    if(p->flagslice4[IJ]>0)
    p->DFBED[IJ]=1;
    
    if(p->flagslice4[IJ]<0)
    p->DFBED[IJ]=-1;
    }
    
    for(int q : dfbed_open)
    p->DFBED[q]=1;
    
    dfbed_inflow(p);
    
    dfbed_comms(p,pgc);
}

void sediment_f::active_zone(lexer *p, ghostcell *pgc)
{
    SLICEBASELOOP
    for(int n=0;n<p->S74;++n)
    if(p->pos_x()>p->S74_xs[n] && p->pos_x()<p->S74_xe[n] && p->pos_y()>p->S74_ys[n] && p->pos_y()<p->S74_ye[n])
    s->active(i,j)=0;

}
void sediment_f::dfbed_inflow(lexer *p)
{
    // SFLOW: no sediment enters through an inflow (B 60, the -i face of a block with G 54): the cells beyond the
    // inflow faces are no sediment cells, which closes the faces for the Exner fluxes, the slope flux and the
    // upwind cells of the lags.  (Before, the face velocity of the sediment was left 0 there, which also halved the
    // velocity and the bed load of the first row, see inflow_PQ.)  The inflow guard S 78 holds the bed of the
    // first row (ini_guard): its bed load is the supply of the cells downstream.  Curvilinear grids only and without
    // a sediment supply at the inflow (S 129 0): the Cartesian grid and S 129 > 0 keep the inflow faces open, the
    // bed load of the ghost cells enters there (fill_PQ_sflow, sediment_exner::qb_inflow).
    if(p->geo_curv!=1 || p->S129!=0)
    return;

    for(int n=0; n<p->gcslin_count; ++n)
    {
        const int ii = p->gcslin[n][0], jj = p->gcslin[n][1], cs = p->gcslin[n][3];
        const int di = (cs==1) ? -1 : (cs==4 ? 1 : 0);
        const int dj = (cs==3) ? -1 : (cs==2 ? 1 : 0);

        for(int q=1; q<=p->margin; ++q)
        {
        const int ia = ii + q*di, ja = jj + q*dj;
            if(ia>=p->imin && ia<p->imin+p->imax && ja>=p->jmin && ja<p->jmin+p->jmax)
            p->DFBED[(ia-p->imin)*p->jmax + (ja-p->jmin)] = -1;
        }
    }
}

void sediment_f::dfbed_comms(lexer *p, ghostcell *pgc)
{
    // sediment cell flag of the neighbour subdomains: the Exner face closure, the bed filter,
    // the non-equilibrium relaxation and the sand slide test DFBED of the ghost cells.
    // Physical boundary ghost cells keep 1 (open), closed walls have zero face velocity.
    pgc->gcslparaxV_int(p,p->DFBED,4);
    pgc->gcslparacoxV_int(p,p->DFBED,1);
}
