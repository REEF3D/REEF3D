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

#include"fnpf_fsf_update.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"slice.h"

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e20)

void fnpf_fsf_update::velcalc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{
    // Per-column version of the former three FLOOP sweeps (bit-identical):
    // grid-metric denominators and the 2D wet-neighbour test are evaluated once
    // per (i,j) column instead of once per cell.
    const int sI = p->jmax*p->kmaxF;
    const int sJ = p->kmaxF;
    const int knoz = p->knoz;
    
    const double *const __restrict Fi = c->Fi;
    double *const __restrict U = c->U;
    double *const __restrict V = c->V;
    double *const __restrict W = c->W;
    
    ILOOP
    JLOOP
    {
        k=0;
        const int c0 = FIJK;
        
        const double dx4 = (-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2]);
        const double dy4 = (-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2]);
        const double dx2 = (p->DXP[IP]+p->DXP[IM1]);
        const double dy2 = (p->DYP[JP]+p->DYP[JM1]);
        const double sz  = p->sigz[IJ];
        
        for(k=0; k<=knoz; ++k)
        {
            const int q = c0 + k;
            
            if(p->flag7[q]>0)
            {
            if(k<knoz)
            {
            U[q] = (-Fi[q+2*sI] + 8.0*Fi[q+sI] - 8.0*Fi[q-sI] + Fi[q-2*sI])/dx4
                 + p->sigx[q]*((Fi[q+1]-Fi[q-1])/(p->DZN[KP]+p->DZN[KM1]));
            
            V[q] = (-Fi[q+2*sJ] + 8.0*Fi[q+sJ] - 8.0*Fi[q-sJ] + Fi[q-2*sJ])/dy4
                 + p->sigy[q]*((Fi[q+1]-Fi[q-1])/(p->DZN[KP]+p->DZN[KM1]));
            
            W[q] = ((-Fi[q+2] + 8.0*Fi[q+1] - 8.0*Fi[q-1] + Fi[q-2])/(-p->ZN[KP2] + 8.0*p->ZN[KP1] - 8.0*p->ZN[KM1] + p->ZN[KM2]))*sz;
            }
            
            if(k==knoz)
            {
            U[q] = (Fi[q+sI]-Fi[q-sI])/dx2
                 + p->sigx[q]*((Fi[q]-Fi[q-1])/(p->DZN[KP]));
            
            V[q] = (Fi[q+sJ]-Fi[q-sJ])/dy2
                 + p->sigy[q]*((Fi[q]-Fi[q-1])/(p->DZN[KP]));
            }
            }
        }
        
        // former FFILOOP4
        if(p->flagslice4[IJ]>0)
        W[c0+knoz] = c->Fz(i,j);
        
        // former third FLOOP: wet test is 2D, clamp only near the inlet
        const bool dryneigh = (p->wet[Im1J]==0 || p->wet[Ip1J]==0 || p->wet[IJm1]==0 || p->wet[IJp1]==0 
                            || p->wet[Im1Jm1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Ip1Jp1]==0
                            || (p->A343>=2 && p->wet[IJ]==0));
        const bool clampcol = (i+p->origin_i<=5);
        
        if(dryneigh || clampcol)
        for(k=0; k<=knoz; ++k)
        {
            const int q = c0 + k;
            
            if(p->flag7[q]>0)
            {
            if(dryneigh)
            {
            U[q]=0.0;
            V[q]=0.0;
            W[q]=0.0;
            }
            
            if(clampcol)
            {
            if(U[q]<=-p->N61)
            U[q] = -0.95*p->N61;
            
            if(U[q]>=p->N61)
            U[q] = 0.95*p->N61;
            
            if(V[q]<=-p->N61)
            V[q] = -0.95*p->N61;
            
            if(V[q]>=p->N61)
            V[q] = 0.95*p->N61;
            
            if(W[q]<=-p->N61)
            W[q] = -0.95*p->N61;
            
            if(W[q]>=p->N61)
            W[q] = 0.95*p->N61;
            }
            }
        }
    }
    
    int gcval=210;
    
    pgc->start7V(p,c->U,c->bc,gcval);
    pgc->start7V(p,c->V,c->bc,gcval);
    pgc->start7V(p,c->W,c->bc,gcval);
    
    
    /*
    k = p->knoz;
    
    SLICELOOP4
    if(fabs(c->W[FIJK])>100.0 && p->YP[JP]>9000 && p->YP[JP]<9400)
    {
    FKLOOP
    {
    cout<<c->Fi[FIJK]<<" "<<((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZP[KP]+p->DZP[KM1]))<<" "<<p->sigz[IJ]<<" "<<c->WL(i,j)<<" "<<p->wd-c->bed(i,j)<<" "<<c->depth(i,j)<<endl;
    }
    cout<<endl;
    }
    
    LOOP
    c->test[IJK] = p->sigz[IJ];
    */
    
    SLICELOOP4
    c->eta_n(i,j) = c->eta(i,j);

    pgc->gcsl_start4(p,c->eta_n,1);   
}

