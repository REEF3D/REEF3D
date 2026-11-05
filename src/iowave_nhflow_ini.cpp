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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::ini_nhflow(lexer *p, fdm_nhf *d, ghostcell* pgc)
{
    wavegen_precalc_ini_nhflow(p,d,pgc);
    wavegen_precalc_relax_func_nhflow(p,pgc);
    
    if(p->B89==1 && p->B98==2)
    nhflow_wavegen_precalc_decomp_space(p,pgc);
    
    if(p->B89==1 && p->B98>=3)
    nhflow_wavegen_precalc_decomp_space_dirichlet(p,pgc);
    
    wavegen_precalc_nhflow(p,d,pgc);
    
    U_relax(p,pgc,d->U,d->UH);
    V_relax(p,pgc,d->V,d->VH);
    W_relax(p,pgc,d->W,d->WH);

    if(p->I30==1)
	full_initialize_nhflow(p,d,pgc);
    
    if(p->B531==1)
    background_initialize_nhflow(p,d,pgc);
}

// B 531 1: the background current of the zones (B 510-514, B 523) in the whole domain at the start, added
// to the velocities of the waves (I 30 1) or of the relaxation zones. A current through a basin at rest
// otherwise has to be spun up from the open edges, which needs a surface slope U/(g t_ramp) over the whole
// basin (1.2 m/s over 8 s: 0.8 m over 70 m). Typical use: a ship in the ship frame, the hull fixed in surge
// and the water passing it at the ship speed (B 514 with t_ramp 0 in B 510). The level of the background
// (eta0) is not set; the background of the first zone that has one is used for the whole domain.
void iowave::background_initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!bg_on)
    {
        if(p->mpirank==0)
        cout<<"B 531: no zone with a background (B 523), the water starts at rest"<<endl;
        return;
    }
    
    // the background of the zones: the generation zones first, then the beaches
    int b=-1;
    
    SLICELOOP4
    {
    if(b<0 && col_gen_bg[IJ]>=0)
    b = col_gen_bg[IJ];
    }
    
    b = pgc->globalimax(b);
    
    if(b<0)
    {
        SLICELOOP4
        {
        if(b<0 && col_beach_bg[IJ]>=0)
        b = col_beach_bg[IJ];
        }
        
        b = pgc->globalimax(b);
    }
    
    if(b<0 || !bgs.carries_current(b))
    {
        if(p->mpirank==0)
        cout<<"B 531: the background of the zones carries no current, the water starts at rest"<<endl;
        return;
    }
    
    double ub0=0.0, vb0=0.0;
    
    LOOP
    WETDRY
    {
        double ub, vb;
        bgs.vel(b,col_h0[IJ],p->XP[IP],p->YP[JP],ub,vb);
        
        if(bgs.profiled(b))
        {
            const double fp = bg_prof(p,b,col_h0[IJ]+bgs.eta(b,p->XP[IP],p->YP[JP]));
            ub *= fp;
            vb *= fp;
        }
        
        d->U[IJK] += ub;
        d->V[IJK] += vb;
        d->UH[IJK] = d->WL(i,j)*d->U[IJK];
        d->VH[IJK] = d->WL(i,j)*d->V[IJK];
        
        ub0 = ub;
        vb0 = vb;
    }
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->UH,10);
    pgc->start4V(p,d->VH,11);
    
    ub0 = pgc->globalmax(fabs(ub0));
    vb0 = pgc->globalmax(fabs(vb0));
    
    if(p->mpirank==0)
    cout<<"B 531: background current in the whole domain at the start, |U| "<<ub0<<" |V| "<<vb0<<" m/s"<<endl;
}