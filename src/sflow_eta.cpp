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

#include"sflow_eta.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"seastate_sflow.h"
#include"hydro_sflow.h"
#include"ioflow.h"
#include"patchBC_interface.h"
#include"sflow_breaking_kb.h"

sflow_eta::sflow_eta(lexer *p, fdm2D *b , ghostcell *pgc, patchBC_interface *ppBC) : eps(1.0e-6), K(p)
{   
    pBC = ppBC;
    
    gcval_eta = 50 + p->F50;     // free surface boundary condition F 50 1-4: ghost-cell labels 51-54
    
    p->phimean=p->F60;
    
    SLICELOOP4
    b->eta(i,j) = 0.0;
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    
    wd_criterion=p->A244;
    
    p->Iarray(temp,p->imax*p->jmax);
    
    if(p->A246==4)
    pkb = new sflow_breaking_kb(p);
}

sflow_eta::~sflow_eta()
{
    delete [] temp;
    delete pkb;
}

void sflow_eta::update(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &WLout, slice &WL0, slice &WLs, double a)
{
    starttime=pgc->timer();
    
    // positivity of the stage with wetting and drying (A 243 > 0) and with REEF3D::HYDRO: the outflow of a cell is at
    // most its water above the film A 244 at the start of the stage plus its rain and sources (donor-cell limiter of
    // the face mass fluxes, conservative). Sheet flow on steep beds and at steps, and the thin layers left on the
    // floodplains after a flood, drained cells below the film: the clip of WL to 0 and the film of wetdry made water
    // (a draining floodplain fed the river with a second discharge for many hours)
    if(b->hydro!=nullptr || p->A243>0)
    positivity_limit(p,b,pgc,WLs);
    
    // continuity: dWL/dt = -div(FE)
    if(p->geo_curv==1)
    {
    SLICELOOP4
    K(i,j) = 0.0;
    
    SLICELOOP4
    WETDRY
    K(i,j) = -(b->FEx(i,j)*p->LEN1[IJ] - b->FEx(i-1,j)*p->LEN1[Im1J]
             + b->FEy(i,j)*p->LEN2[IJ] - b->FEy(i,j-1)*p->LEN2[IJm1])/p->AREA2D[IJ];
    }
    
    else
    {
    const int nj = p->jmax;
    const int *const flag = p->flagslice4, *const wet = p->wet;
    const double *const FEx = b->FEx.V, *const FEy = b->FEy.V;
    double *const Kv = K.V;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    const double dxn = p->DXN[ii+marge];
    
        for(int jj=0; jj<p->knoy; ++jj)
        {
        const int q = q0 + jj;
        
        if(flag[q]>0)
        Kv[q] = wet[q]==1 ? -(FEx[q] - FEx[q-nj])/dxn
                            -(FEy[q] - FEy[q-1])/p->DYN[jj+marge]*p->y_dir
                          : 0.0;
        }
    }
    }

    // REEF3D::SEASTATE vortex force (A 751 2): Stokes transport
    if(b->wave!=nullptr)
    b->wave->mass_source(p,b,K);

    // REEF3D::HYDRO: rain and point sources, on wet and dry cells
    if(b->hydro!=nullptr)
    b->hydro->mass_source(p,b,K);
    
    // depth-averaged suspended load (S 120 > 0): the mass flux of the step, weighted as the stages, so that
    // WL^n+1 = WL^n - dt div(FEacc) (stage output a WL0 + (1-a)(WLs + dt K): FEacc = (1-a)(FEacc + FE))
    if(b->FEacc_x!=nullptr && p->geo_curv==0)
    {
    const double w = 1.0-a;
    
        if(a==0.0)
        {
        for(i=-1; i<p->knox; ++i)
        for(j=0; j<p->knoy; ++j)
        (*b->FEacc_x)(i,j) = 0.0;
        
        for(i=0; i<p->knox; ++i)
        for(j=-1; j<p->knoy; ++j)
        (*b->FEacc_y)(i,j) = 0.0;
        }
    
        for(i=-1; i<p->knox; ++i)
        for(j=0; j<p->knoy; ++j)
        (*b->FEacc_x)(i,j) = w*((*b->FEacc_x)(i,j) + b->FEx(i,j));
        
        for(i=0; i<p->knox; ++i)
        for(j=-1; j<p->knoy; ++j)
        (*b->FEacc_y)(i,j) = w*((*b->FEacc_y)(i,j) + b->FEy(i,j));
    
    b->FEacc_step = p->count;
    }
    
    // eta + boundary conditions
    {
    const int *const flag = p->flagslice4;
    const double *const w0 = WL0.V, *const ws = WLs.V, *const Kv = K.V, *const depth = b->depth.V;
    double *const wo = WLout.V, *const eta = b->eta.V;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*p->jmax - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag[q]>0)
        {
        wo[q] = a*w0[q] + (1.0-a)*(ws[q] + p->dt*Kv[q]);
        eta[q] = wo[q] - depth[q];
        }
    }
    }
    
    
    pflow->waterlevel2D(p,b,pgc,b->eta);
    pflow->eta_relax(p,pgc,b->eta);
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    
    {
    const int *const flag = p->flagslice4;
    const double *const depth = b->depth.V, *const eta = b->eta.V;
    double *const wo = WLout.V;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*p->jmax - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag[q]>0)
        wo[q] = MAX(eta[q] + depth[q], 0.0);
    }
    }
    
    pgc->gcsl_start4(p,WLout,gcval_eta);
    
    // wetting & drying
    wetdry(p,b,pgc,WLout);
    
    p->lsmtime+=pgc->timer()-starttime;
}

void sflow_eta::positivity_limit(lexer *p, fdm2D *b, ghostcell *pgc, slice &WLs)
{
    const double dt = p->dt;
    double out;

    // the factor of the outflow of each cell (in K, which update sets afterwards)
    // Cartesian grid without REEF3D::HYDRO: over the arrays of the slices
    if(p->geo_curv!=1 && b->hydro==nullptr)
    {
    const int nj = p->jmax;
    const int *const flag = p->flagslice4;
    const double *const FEx = b->FEx.V, *const FEy = b->FEy.V, *const wls = WLs.V;
    double *const Kv = K.V;
    const double wdc = wd_criterion;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    const double dxn = p->DXN[ii+marge];
    
        for(int jj=0; jj<p->knoy; ++jj)
        {
        const int q = q0 + jj;
        
        if(flag[q]>0)
        {
        const double fe = FEx[q], fw = FEx[q-nj], fn = FEy[q], fs = FEy[q-1];
        
        const double o = (MAX(fe,0.0) + MAX(-fw,0.0))/dxn
                       + (MAX(fn,0.0) + MAX(-fs,0.0))/p->DYN[jj+marge]*p->y_dir;
        
        // the film A 244 stays: wetdry gives a cell that drains below it the film again (water from nothing)
        const double avail = MAX(wls[q]-wdc,0.0) + 0.0;
        
        Kv[q] = (dt*o>avail && o>0.0) ? avail/(dt*o) : 1.0;
        }
        }
    }
    }
    
    else
    SLICELOOP4
    {
        if(p->geo_curv==1)
        out = (MAX(b->FEx(i,j),0.0)*p->LEN1[IJ] + MAX(-b->FEx(i-1,j),0.0)*p->LEN1[Im1J]
            +  MAX(b->FEy(i,j),0.0)*p->LEN2[IJ] + MAX(-b->FEy(i,j-1),0.0)*p->LEN2[IJm1])/p->AREA2D[IJ];
        else
        out = (MAX(b->FEx(i,j),0.0) + MAX(-b->FEx(i-1,j),0.0))/p->DXN[IP]
            + (MAX(b->FEy(i,j),0.0) + MAX(-b->FEy(i,j-1),0.0))/p->DYN[JP]*p->y_dir;
        // the film A 244 stays: wetdry gives a cell that drains below it the film again (water from nothing)
        const double avail = MAX(WLs(i,j)-wd_criterion,0.0) + (b->hydro!=nullptr ? dt*b->hydro->S(i,j) : 0.0);
        K(i,j) = (dt*out>avail && out>0.0) ? avail/(dt*out) : 1.0;
    }

    // the neighbours across the ranks; the ghost cells of the domain boundaries do not limit (inflow)
    pgc->gcsl_start4(p,K,50);

    for(n=0;n<p->gcbsl4_count;++n)
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];
    const int cs = p->gcbsl4[n][3];

    K(i-(cs==1?1:0)+(cs==4?1:0), j-(cs==3?1:0)+(cs==2?1:0)) = 1.0;
    }

    // a face takes the factor of its donor cell
    {
    const int nj = p->jmax;
    const int *const flag1 = p->flagslice1, *const flag2 = p->flagslice2;
    const double *const Kv = K.V;
    double *const FEx = b->FEx.V, *const FEy = b->FEy.V;
    const int ie = p->knox - p->ulast;
    const int je = p->knoy - p->vlast;
    
    for(int ii=0; ii<ie; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag1[q]>0)
        FEx[q] *= FEx[q]>0.0 ? Kv[q] : Kv[q+nj];
    }
    
    if(p->j_dir==1)
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+je; ++q)
        if(flag2[q]>0)
        FEy[q] *= FEy[q]>0.0 ? Kv[q] : Kv[q+1];
    }
    }
}

void sflow_eta::depth_update(lexer *p, fdm2D *b , ghostcell *pgc, slice &WL)
{
	// still water depth
	SLICELOOP4
	b->depth(i,j) = p->wd - b->bed(i,j);
    
    pgc->gcsl_start4(p,b->depth,50);
    
    // the water column is kept: eta follows the bed
    SLICELOOP4
    b->eta(i,j) = WL(i,j) - b->depth(i,j);
    
    // set fsf outflow
    double wsfout=p->phimean;
    double f=1.0;
    
    // F 62 given (unset -1.0e20): also a level at or below the datum z = 0
    if(p->F62>-1.0e20)
    {
        if(p->F64==0)
        wsfout=p->F62;
        
        if(p->F64>0)
        {
        if(p->count<p->F64)
        f = 0.5*cos(PI + PI*double(p->count)/double(p->F64)) + 0.5;
        
        if(p->count>=p->F64)
        f = 1.0;
        
        wsfout = f*p->F62 + (1.0-f)*p->F60;
        }
    }
    
    if(p->F50==2 || p->F50==3)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
    // not at iowave Riemann / Flather / clamped edges (B 520 methods 3-6)
    const int cs = p->gcslout[n][3];
    if((cs==4 && p->open_xp==1) || (cs==3 && p->open_ym==1) || (cs==2 && p->open_yp==1) || (cs==1 && p->open_xm==1))
    continue;
    
    if(p->wet[IJ]==1)
    b->eta(i,j) = wsfout-p->wd;
    }
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    
    SLICELOOP4
    WL(i,j) = MAX(b->eta(i,j) + b->depth(i,j), 0.0);
    
    pgc->gcsl_start4(p,WL,gcval_eta);
    
    // wetdry
    wetdry(p,b,pgc,WL);
}

void sflow_eta::ini(lexer *p, fdm2D *b , ghostcell *pgc, ioflow *pflow)
{
    p->phimean=p->F60;
    
    SLICELOOP4
    b->eta(i,j) = 0.0;
    
    pflow->eta_relax(p,pgc,b->eta);
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
}

void sflow_eta::breaking_zones(lexer *p, ioflow *pflow)
{
    if(pkb!=nullptr)
    pkb->zones(p,pflow);
}
