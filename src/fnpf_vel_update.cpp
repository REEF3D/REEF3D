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
#include"fnpf_sigma_col.h"

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
    // surface only (no 3D U, V, W, fnpf_memory.h): the node k = knoz of every column, as in the
    // loops below, into Us, Vs, Ws; their halo as start7V(210) makes it at that node: the partition
    // (with the corners), then the wall copies of fivec_vel
    if(c->U==nullptr)
    {
        ILOOP
        JLOOP
        {
            k=0;
            const int c0 = FIJK;
            const int q = c0 + knoz;
            
            k=knoz;
            
            if(p->flag7[q]>0)
            {
            const double dx2 = (p->DXP[IP]+p->DXP[IM1]);
            const double dy2 = (p->DYP[JP]+p->DYP[JM1]);
            
            fnpf_sigma_col col;
            
            if(p->sigx==nullptr)
            col.set(p,c,i,j);
            
            const double sgx = p->sigx!=nullptr ? p->sigx[q] : col.sx(p->ZN[KP]);
            const double sgy = p->sigy!=nullptr ? p->sigy[q] : col.sy(p->ZN[KP]);
            
            c->Us(i,j) = (Fi[q+sI]-Fi[q-sI])/dx2
                       + sgx*((Fi[q]-Fi[q-1])/(p->DZN[KP]));
            
            c->Vs(i,j) = (Fi[q+sJ]-Fi[q-sJ])/dy2
                       + sgy*((Fi[q]-Fi[q-1])/(p->DZN[KP]));
            }
            
            if(p->flagslice4[IJ]>0)
            c->Ws(i,j) = c->Fz(i,j);
            
            const bool dryneigh = (p->wet[Im1J]==0 || p->wet[Ip1J]==0 || p->wet[IJm1]==0 || p->wet[IJp1]==0 
                                || p->wet[Im1Jm1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Ip1Jp1]==0
                                || (p->A343>=2 && p->wet[IJ]==0));
            
            const bool clampcol = (i+p->origin_i<=5);
            
            if((dryneigh || clampcol) && p->flag7[q]>0)
            {
                if(dryneigh)
                {
                c->Us(i,j)=0.0;
                c->Vs(i,j)=0.0;
                c->Ws(i,j)=0.0;
                }
                
                if(clampcol)
                {
                if(c->Us(i,j)<=-p->N61)
                c->Us(i,j) = -0.95*p->N61;
                if(c->Us(i,j)>=p->N61)
                c->Us(i,j) = 0.95*p->N61;
                if(c->Vs(i,j)<=-p->N61)
                c->Vs(i,j) = -0.95*p->N61;
                if(c->Vs(i,j)>=p->N61)
                c->Vs(i,j) = 0.95*p->N61;
                if(c->Ws(i,j)<=-p->N61)
                c->Ws(i,j) = -0.95*p->N61;
                if(c->Ws(i,j)>=p->N61)
                c->Ws(i,j) = 0.95*p->N61;
                }
            }
        }
        
        slice *const vs[3] = {&c->Us,&c->Vs,&c->Ws};
        
        for(slice *v : vs)
        {
            pgc->gcslparax(p,*v,4);
            pgc->gcslparacox(p,*v,4);
            
            k=knoz;
            ILOOP
            JLOOP
            if(p->flag7[FIJK]>0)
            {
                slice &g = *v;
                
                if(p->flag7[FIm1JK]<0)
                g(i-1,j) = g(i-2,j) = g(i-3,j) = g(i,j);
                
                if(p->flag7[FIp1JK]<0)
                g(i+1,j) = g(i+2,j) = g(i+3,j) = g(i,j);
                
                if(p->flag7[FIJm1K]<0)
                g(i,j-1) = g(i,j-2) = g(i,j-3) = g(i,j);
                
                if(p->flag7[FIJp1K]<0)
                g(i,j+1) = g(i,j+2) = g(i,j+3) = g(i,j);
            }
        }
        
        c->vel_state = 1;
        
        SLICELOOP4
        c->eta_n(i,j) = c->eta(i,j);
        
        pgc->gcsl_start4(p,c->eta_n,1);
        
        return;
    }
    
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
        
        // the sigma metrics of the column (without the 3D arrays, fnpf_sigma_col.h)
        fnpf_sigma_col col;
        
        if(p->sigx==nullptr)
        col.set(p,c,i,j);
        
        for(k=0; k<=knoz; ++k)
        {
            const int q = c0 + k;
            
            if(p->flag7[q]>0)
            {
            const double sgx = p->sigx!=nullptr ? p->sigx[q] : col.sx(p->ZN[KP]);
            const double sgy = p->sigy!=nullptr ? p->sigy[q] : col.sy(p->ZN[KP]);
            
            if(k<knoz)
            {
            U[q] = (-Fi[q+2*sI] + 8.0*Fi[q+sI] - 8.0*Fi[q-sI] + Fi[q-2*sI])/dx4
                 + sgx*((Fi[q+1]-Fi[q-1])/(p->DZN[KP]+p->DZN[KM1]));
            
            V[q] = (-Fi[q+2*sJ] + 8.0*Fi[q+sJ] - 8.0*Fi[q-sJ] + Fi[q-2*sJ])/dy4
                 + sgy*((Fi[q+1]-Fi[q-1])/(p->DZN[KP]+p->DZN[KM1]));
            
            W[q] = ((-Fi[q+2] + 8.0*Fi[q+1] - 8.0*Fi[q-1] + Fi[q-2])/(-p->ZN[KP2] + 8.0*p->ZN[KP1] - 8.0*p->ZN[KM1] + p->ZN[KM2]))*sz;
            }
            
            if(k==knoz)
            {
            U[q] = (Fi[q+sI]-Fi[q-sI])/dx2
                 + sgx*((Fi[q]-Fi[q-1])/(p->DZN[KP]));
            
            V[q] = (Fi[q+sJ]-Fi[q-sJ])/dy2
                 + sgy*((Fi[q]-Fi[q-1])/(p->DZN[KP]));
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
    
    c->vel_state = 1;
    
    int gcval=210;
    
    pgc->start7V(p,c->U,c->bc,gcval);
    pgc->start7V(p,c->V,c->bc,gcval);
    pgc->start7V(p,c->W,c->bc,gcval);
    
    
    SLICELOOP4
    c->eta_n(i,j) = c->eta(i,j);

    pgc->gcsl_start4(p,c->eta_n,1);   
}

