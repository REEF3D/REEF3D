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

#include"sflow_HLL.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"sflow_flux_build_f.h"

sflow_HLL::sflow_HLL(lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new sflow_flux_build_f(p,pgc,pBC);
    
    // physical in- and outflow faces carry the flux of the ghost state,
    // all other boundary faces are walls
    inflow=0;
    outflow=0;
    
    if(p->B98>=3 || p->B60>=1)
    inflow=1;
    
    if(p->B99>=3 || p->B60>=1)
    outflow=1;
}

sflow_HLL::~sflow_HLL()
{
}

void sflow_HLL::start(lexer *p, fdm2D *b, int ipol)
{
    if(ipol==1)
    aij_U(p,b);

    if(ipol==2 && p->j_dir==1)
    aij_V(p,b);

    if(ipol==3)
    aij_W(p,b);
    
    if(ipol==4)
    aij_E(p,b);
}

void sflow_HLL::aij_U(lexer *p, fdm2D *b)
{
    pflux->start_U(p,b,pgc);
    
    // Boussinesq: fluxes from the volume flux M, dissipation from the conserved V
    if(p->A220==4)
    HLL(p,b,b->QUs,b->QUn,b->QUe,b->QUw,b->Fx,b->Fy);
    
    else
    HLL(p,b,b->UHs,b->UHn,b->UHe,b->UHw,b->Fx,b->Fy);
    flux_bc(p,b,1);
    
    divergence(p,b,b->F);
}

void sflow_HLL::aij_V(lexer *p, fdm2D *b)
{
    pflux->start_V(p,b,pgc);
    
    if(p->A220==4)
    HLL(p,b,b->QVs,b->QVn,b->QVe,b->QVw,b->Fx,b->Fy);
    
    else
    HLL(p,b,b->VHs,b->VHn,b->VHe,b->VHw,b->Fx,b->Fy);
    flux_bc(p,b,2);
    
    divergence(p,b,b->G);
}

void sflow_HLL::aij_W(lexer *p, fdm2D *b)
{
    pflux->start_W(p,b,pgc);
    HLL(p,b,b->WHs,b->WHn,b->WHe,b->WHw,b->Fx,b->Fy);
    flux_bc(p,b,3);
    
    divergence(p,b,b->H);
}

void sflow_HLL::aij_E(lexer *p, fdm2D *b)
{
    pflux->start_E(p,b,pgc);
    HLL(p,b,b->Ds,b->Dn,b->De,b->Dw,b->FEx,b->FEy);
    
    // no mass flux into or out of dry cells
    SLICELOOP4
    WETDRY
    {
    if(p->wet[Ip1J]==0)
    b->FEx(i,j) = 0.0;
    
    if(p->wet[Im1J]==0)
    b->FEx(i-1,j) = 0.0;
    
    if(p->wet[IJp1]==0)
    b->FEy(i,j) = 0.0;
    
    if(p->wet[IJm1]==0)
    b->FEy(i,j-1) = 0.0;
    }
    
    flux_bc(p,b,4);
}

void sflow_HLL::divergence(lexer *p, fdm2D *b, slice &f)
{
    SLICELOOP4
    WETDRY
    f(i,j) -= (b->Fx(i,j) - b->Fx(i-1,j))/p->DXN[IP] 
            + (b->Fy(i,j) - b->Fy(i,j-1))/p->DYN[JP]*p->y_dir;
}

void sflow_HLL::HLL(lexer *p, fdm2D *b, slice &Qs, slice &Qn, slice &Qe, slice &Qw, slice &Fx, slice &Fy)
{    
    // HLL flux x-dir
    SLICELOOP1
    {
        if(b->Ss(i,j)>=0.0)
        Fx(i,j) = b->Fs(i,j);
        
        else
        if(b->Sn(i,j)<=0.0)
        Fx(i,j) = b->Fn(i,j);
        
        else
        {
        denom = b->Sn(i,j)-b->Ss(i,j);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fx(i,j) = (b->Sn(i,j)*b->Fs(i,j) - b->Ss(i,j)*b->Fn(i,j) + b->Sn(i,j)*b->Ss(i,j)*(Qn(i,j) - Qs(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    SLICELOOP2
    {
        if(b->Se(i,j)>=0.0)
        Fy(i,j) = b->Fe(i,j);
        
        else
        if(b->Sw(i,j)<=0.0)
        Fy(i,j) = b->Fw(i,j);
        
        else
        {
        denom = b->Sw(i,j)-b->Se(i,j);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fy(i,j) = (b->Sw(i,j)*b->Fe(i,j) - b->Se(i,j)*b->Fw(i,j) + b->Sw(i,j)*b->Se(i,j)*(Qw(i,j) - Qe(i,j)))/denom;
        }
    }
}

void sflow_HLL::flux_bc(lexer *p, fdm2D *b, int ipol)
{
    // MPI: the first face of each partition is computed by the neighbour
    // (exchange only, the physical boundary faces are set below)
    slice &Fxp = (ipol==4)?b->FEx:b->Fx;
    slice &Fyp = (ipol==4)?b->FEy:b->Fy;
    
    if(p->mpi_size>1)
    {
    pgc->gcslparax(p,Fxp,1);
    pgc->gcslparacox(p,Fxp,14);
    
        if(p->j_dir==1)
        {
        pgc->gcslparax(p,Fyp,2);
        pgc->gcslparacox(p,Fyp,14);
        }
    }
    
    // boundary faces
    //  walls: no mass flux, momentum flux = hydrostatic pressure of the wall cell
    //  inflow (gcslin) / outflow (gcslout): physical flux of the ghost state
    slice &Fx = (ipol==4)?b->FEx:b->Fx;
    slice &Fy = (ipol==4)?b->FEy:b->Fy;
    
    const double g = fabs(p->W22);
    
    SLICELOOP4
    {
        // x-dir
        if(p->flagslice4[Im1J]<0)
        {
        if(ipol==1)
        Fx(i-1,j) = 0.5*g*b->eta(i,j)*b->eta(i,j) + g*b->eta(i,j)*b->dfx(i-1,j);
        
        else
        Fx(i-1,j) = 0.0;
        }
        
        if(p->flagslice4[Ip1J]<0)
        {
        if(ipol==1)
        Fx(i,j) = 0.5*g*b->eta(i,j)*b->eta(i,j) + g*b->eta(i,j)*b->dfx(i,j);
        
        else
        Fx(i,j) = 0.0;
        }
        
        // y-dir
        if(p->j_dir==1)
        {
            if(p->flagslice4[IJm1]<0)
            {
            if(ipol==2)
            Fy(i,j-1) = 0.5*g*b->eta(i,j)*b->eta(i,j) + g*b->eta(i,j)*b->dfy(i,j-1);
            
            else
            Fy(i,j-1) = 0.0;
            }
            
            if(p->flagslice4[IJp1]<0)
            {
            if(ipol==2)
            Fy(i,j) = 0.5*g*b->eta(i,j)*b->eta(i,j) + g*b->eta(i,j)*b->dfy(i,j);
            
            else
            Fy(i,j) = 0.0;
            }
        }
    }
    
    // inflow: ghost state from the ghost velocity and water level
    double wl;
    
    if(inflow==1)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0]-1;
    j=p->gcslin[n][1];
    
    wl = MAX(b->eta(i,j) + b->depth(i,j), 0.0);
    
        if(ipol==1)
        Fx(i,j) = wl*b->U(i,j)*b->U(i,j) + 0.5*g*b->eta(i,j)*b->eta(i,j) + g*b->eta(i,j)*b->dfx(i,j);
        
        if(ipol==2)
        Fx(i,j) = wl*b->V(i,j)*b->U(i,j);
        
        if(ipol==3)
        Fx(i,j) = wl*b->W(i,j)*b->U(i,j);
        
        if(ipol==4)
        Fx(i,j) = wl*b->U(i,j);
    }
    
    // patch boundaries: physical flux of the ghost state
    int cs;
    double eg,ug,vg,wg;
    
    GCSL4LOOP
    if(p->gcbsl4[n][4]>=100)
    {
    i  = p->gcbsl4[n][0];
    j  = p->gcbsl4[n][1];
    cs = p->gcbsl4[n][3];
    
    int ii=0,jj=0;
    
    if(cs==1)
    ii=-1;
    
    if(cs==4)
    ii=1;
    
    if(cs==2)
    jj=1;
    
    if(cs==3)
    jj=-1;
    
    eg = b->eta(i+ii,j+jj);
    wl = MAX(eg + b->depth(i+ii,j+jj), 0.0);
    ug = b->U(i+ii,j+jj);
    vg = b->V(i+ii,j+jj);
    wg = b->W(i+ii,j+jj);
    
        // x-face
        if(cs==1 || cs==4)
        {
        int fi = (cs==1)?i-1:i;
        
        if(ipol==1)
        Fx(fi,j) = wl*ug*ug + 0.5*g*eg*eg + g*eg*b->dfx(fi,j);
        
        if(ipol==2)
        Fx(fi,j) = wl*vg*ug;
        
        if(ipol==3)
        Fx(fi,j) = wl*wg*ug;
        
        if(ipol==4)
        Fx(fi,j) = wl*ug;
        }
        
        // y-face
        if((cs==2 || cs==3) && p->j_dir==1)
        {
        int fj = (cs==3)?j-1:j;
        
        if(ipol==1)
        Fy(i,fj) = wl*ug*vg;
        
        if(ipol==2)
        Fy(i,fj) = wl*vg*vg + 0.5*g*eg*eg + g*eg*b->dfy(i,fj);
        
        if(ipol==3)
        Fy(i,fj) = wl*wg*vg;
        
        if(ipol==4)
        Fy(i,fj) = wl*vg;
        }
    }
    
    // outflow
    if(outflow==1)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
    wl = MAX(b->eta(i+1,j) + b->depth(i+1,j), 0.0);
    
        if(ipol==1)
        Fx(i,j) = wl*b->U(i+1,j)*b->U(i+1,j) + 0.5*g*b->eta(i+1,j)*b->eta(i+1,j) + g*b->eta(i+1,j)*b->dfx(i,j);
        
        if(ipol==2)
        Fx(i,j) = wl*b->V(i+1,j)*b->U(i+1,j);
        
        if(ipol==3)
        Fx(i,j) = wl*b->W(i+1,j)*b->U(i+1,j);
        
        if(ipol==4)
        Fx(i,j) = wl*b->U(i+1,j);
    }
}
