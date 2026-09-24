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

#include"sflow_momentum_func.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"sflow_fsf.h"
#include"sflow_signal_speed.h"
#include"sflow_reconstruct.h"
#include"sflow_HLL.h"
#include"sflow_diffusion.h"
#include"sflow_pressure.h"
#include"sflow_forcing.h"
#include"sflow_rough_manning.h"
#include"sflow_rough_void.h"
#include"sflow_rheology_f.h"
#include"sflow_rheology_v.h"
#include"solver2D.h"
#include"6DOF.h"

#define WLVL (fabs(WL(i,j))>p->A244?WL(i,j):1.0e20)

sflow_momentum_func::sflow_momentum_func(lexer *p, fdm2D *b, ghostcell *pgc, sflow_HLL *pphll, sflow_signal_speed *ppss, 
                                         sflow_reconstruct *pprecon, sflow_diffusion *ppdiff, sflow_pressure *pppress, 
                                         solver2D *ppsolv, solver2D *pppoissonsolv, ioflow *ppflow, sflow_fsf *ppfsf, 
                                         sflow_forcing *ppsfdf, sixdof *pp6dof)
                                        : UHDIFF(p),VHDIFF(p),WHDIFF(p),Un(p),Vn(p),etaS(p)
{
    phll = pphll;
    pss = ppss;
    precon = pprecon;
    pdiff = ppdiff;
    ppress = pppress;
    psolv = ppsolv;
    ppoissonsolv = pppoissonsolv;
    pflow = ppflow;
    pfsf = ppfsf;
    psfdf = ppsfdf;
    p6dof = pp6dof;
    
    if(p->A218==0)
    prough = new sflow_rough_void(p);
    
    if(p->A218==1)
    prough = new sflow_rough_manning(p);
    
    if(p->W90==0)
    prheo = new sflow_rheology_v(p);
    
    if(p->W90==1)
    prheo = new sflow_rheology_f(p);
    

	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    // conserved variables: MPI exchange only, physical ghost cells are set in ghostcells()
    gcval_uh=10;
	gcval_vh=10;
	gcval_wh=10;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    inflow_flag=0;
    outflow_flag=0;
    
    if(p->B98>=3 || p->B60>=1)
    inflow_flag=1;
    
    if(p->B99>=3)
    outflow_flag=1;
    
    if(p->B60>=1)
    outflow_flag=2;
}

sflow_momentum_func::~sflow_momentum_func()
{
}

void sflow_momentum_func::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // conserved variables from the initial velocities and water depth
    SLICELOOP4
    {
    if(p->wet[IJ]==1)
    {
    b->UH(i,j) = b->U(i,j)*b->WL(i,j);
    b->VH(i,j) = b->V(i,j)*b->WL(i,j)*p->y_dir;
    b->WH(i,j) = b->W(i,j)*b->WL(i,j);
    }
    
    if(p->wet[IJ]==0)
    {
    b->UH(i,j) = 0.0;
    b->VH(i,j) = 0.0;
    b->WH(i,j) = 0.0;
    }
    }
    
    inflow(p,b,pgc,pflow);
    
    velcalc(p,b,pgc,b->UH,b->VH,b->WH,b->WL,2);
}

void sflow_momentum_func::inflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
    // ghost cell velocities at in- and outflow boundaries (cell centred)
    pflow->discharge2D(p,b,pgc);
    pflow->inflow2D(p,b,pgc,b->U,b->V,b->bed,b->eta);
    
    // Dirichlet wave generation also provides w
    if(p->B98>=3)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        for(q=1;q<=3;++q)
        b->W(i-q,j) = b->ws(i-q,j);
    }
    
    if(p->B60>=1 && p->B98<3)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        for(q=1;q<=3;++q)
        b->W(i-q,j) = 0.0;
    }
    
    ghostcells(p,b,pgc,b->UH,b->VH,b->WH,b->WL);
}

void sflow_momentum_func::stage(lexer *p, fdm2D *b, ghostcell *pgc, 
                                slice &WLs, slice &UHs, slice &VHs, slice &WHs, 
                                slice &WLo, slice &UHo, slice &VHo, slice &WHo, double a, int iter, bool finalize)
{
    const double alpha = 1.0-a;
    
    // stage state; Un,Vn: velocity the RK update starts from, a*u_n + (1-a)*u_s
    // (used for Du/Dt in the quadratic pressure)
    double wln;
    
    SLICELOOP4
    {
    wln = b->WL(i,j)>p->A244?b->WL(i,j):1.0e20;
    
    Un(i,j) = a*b->UH(i,j)/wln + alpha*b->U(i,j);
    Vn(i,j) = a*b->VH(i,j)/wln + alpha*b->V(i,j);
    etaS(i,j) = b->eta(i,j);
    }
    
    reconstruct(p,b,pgc,WLs,UHs,VHs,WHs);
    
    // FSF flux
    starttime=pgc->timer();
    phll->start(p,b,4);
    p->fsftime+=pgc->timer()-starttime;
    
    clearrhs(p,b);
    
    // U
	starttime=pgc->timer();
    phll->start(p,b,1);
	ppress->upgrad(p,b,b->eta);
    prough->u_source(p,b,b->U,b->V,WLs);
    prheo->u_source(p,b,b->U,b->V,WLs);
    p6dof->isource2D(p,b,pgc);
	irhs(p,b);
	pdiff->diff_u(p,b,pgc,psolv,UHDIFF,UHs,b->U,b->V,WLs,alpha);
    p->utime+=pgc->timer()-starttime;
    
	// V
	starttime=pgc->timer();
    if(p->j_dir==1)
    {
    phll->start(p,b,2);
	ppress->vpgrad(p,b,b->eta);
    prough->v_source(p,b,b->U,b->V,WLs);
    prheo->v_source(p,b,b->U,b->V,WLs);
    p6dof->jsource2D(p,b,pgc);
	jrhs(p,b);
    }
	pdiff->diff_v(p,b,pgc,psolv,VHDIFF,VHs,b->U,b->V,WLs,alpha);
    p->vtime+=pgc->timer()-starttime;
    
    // W
    if(p->A220>0)
    {
    if(p->A214==1)
    phll->start(p,b,3);
    ppress->wpgrad(p,b,b->eta);
    krhs(p,b);
    pdiff->diff_w(p,b,pgc,psolv,WHDIFF,WHs,b->U,b->V,b->W,WLs,alpha);
    }
    
    // update FSF: WL, eta, wetdry
    pfsf->update(p,b,pgc,pflow,WLo,b->WL,WLs,a);
    pfsf->breaking(p,b,pgc,b->eta,etaS,alpha);
    
    // update momentum
    SLICELOOP4
    {
	UHo(i,j) = a*b->UH(i,j) + alpha*(UHDIFF(i,j) + p->dt*b->F(i,j));
    VHo(i,j) = (a*b->VH(i,j) + alpha*(VHDIFF(i,j) + p->dt*b->G(i,j)))*p->y_dir;
    }
    
    if(p->A220>0)
    SLICELOOP4
    WHo(i,j) = a*b->WH(i,j) + alpha*(WHDIFF(i,j) + p->dt*b->H(i,j));
    
    if(p->A220==0)
    SLICELOOP4
    WHo(i,j) = 0.0;
    
    velcalc(p,b,pgc,UHo,VHo,WHo,WLo,1);
    
    // direct forcing
    psfdf->forcing(p,b,pgc,p6dof,iter,alpha,UHo,VHo,WHo,WLo,finalize);
    
    // non-hydrostatic pressure
    ppress->start(p,b,pgc,ppoissonsolv,pflow,UHo,VHo,WHo,WLo,Un,Vn,alpha);
    
    velcalc(p,b,pgc,UHo,VHo,WHo,WLo,0);
    
    // relaxation zones
    pflow->um_relax(p,pgc,b->U,UHo,WLo);
    pflow->vm_relax(p,pgc,b->V,VHo,WLo);
    pflow->wm_relax(p,pgc,b->W,WHo,WLo);
    pflow->pm_relax(p,pgc,b->press);
    
    velcalc(p,b,pgc,UHo,VHo,WHo,WLo,2);
}

void sflow_momentum_func::reconstruct(lexer *p, fdm2D *b, ghostcell *pgc, slice &WL, slice &UH, slice &VH, slice &WH)
{
    starttime=pgc->timer();
    
    // eta and face water depth
    precon->reconstruct_x(p, pgc, b, b->eta, b->ETAs, b->ETAn);
    precon->reconstruct_y(p, pgc, b, b->eta, b->ETAe, b->ETAw);
    precon->reconstruct_WL(p, pgc, b);
    
    // velocities: only the normal components enter the signal speeds and fluxes
    precon->reconstruct_x(p, pgc, b, b->U, b->Us, b->Un);
    precon->reconstruct_y(p, pgc, b, b->V, b->Ve, b->Vw);

    // conserved variables
    precon->reconstruct_x(p, pgc, b, UH, b->UHs, b->UHn);
    precon->reconstruct_y(p, pgc, b, UH, b->UHe, b->UHw);

    if(p->j_dir==1)
    {
    precon->reconstruct_x(p, pgc, b, VH, b->VHs, b->VHn);
    precon->reconstruct_y(p, pgc, b, VH, b->VHe, b->VHw);
    }

    if(p->A220>0 && p->A214==1)
    {
    precon->reconstruct_x(p, pgc, b, WH, b->WHs, b->WHn);
    precon->reconstruct_y(p, pgc, b, WH, b->WHe, b->WHw);
    }
    
    // wetdry
    pfsf->wetdry_fluxes(p,b,pgc,WL);
    
    // signal speeds
    pss->signal_speed_update(p, pgc, b, b->Us, b->Un, b->Ve, b->Vw, b->Ds, b->Dn, b->De, b->Dw);
    
    p->recontime+=pgc->timer()-starttime;
}

void sflow_momentum_func::velcalc(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WH, slice &WL, int mode)
{
    // mode 0: interior cells, 1: + ghost cells, 2: + ghost cells + staggered diagnostics
    const double g = fabs(p->W22);
    double lim;
    
    SLICELOOP4
    {
        // Froude number limiter
        if(p->wet[IJ]==1)
        {
        lim = p->A531*WL(i,j)*sqrt(g*WL(i,j));
        
        UH(i,j) = MAX(MIN(UH(i,j), lim), -lim);
        VH(i,j) = MAX(MIN(VH(i,j), lim), -lim);
        WH(i,j) = MAX(MIN(WH(i,j), lim), -lim);
        }
        
        if(p->B60>=1)
        if(p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0)
        {
        lim = 0.1*p->A531*WL(i,j)*sqrt(g*WL(i,j));
        
        UH(i,j) = MAX(MIN(UH(i,j), lim), -lim);
        VH(i,j) = MAX(MIN(VH(i,j), lim), -lim);
        WH(i,j) = MAX(MIN(WH(i,j), lim), -lim);
        }
        
        // velocities
        if(p->wet[IJ]==1)
        {
        b->U(i,j) = UH(i,j)/WLVL;
        b->V(i,j) = VH(i,j)/WLVL*p->y_dir;
        b->W(i,j) = WH(i,j)/WLVL;
        }
        
        if(p->wet[IJ]==0)
        {
        UH(i,j) = 0.0;
        VH(i,j) = 0.0;
        WH(i,j) = 0.0;
        
        b->U(i,j) = 0.0;
        b->V(i,j) = 0.0;
        b->W(i,j) = 0.0;
        }
    }
    
    if(mode>=1)
    ghostcells(p,b,pgc,UH,VH,WH,WL);
    
    if(mode>=2)
    face_velocities(p,b,pgc,WL);
}

void sflow_momentum_func::mpi4(lexer *p, ghostcell *pgc, slice &f)
{
    // MPI exchange only, physical ghost cells are set in ghostcells()
    if(p->mpi_size>1)
    {
    pgc->gcslparax(p,f,4);
    pgc->gcslparacox(p,f,gcval_u);
    }
}

void sflow_momentum_func::ghostcells(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WH, slice &WL)
{
    int cs,bc;
    double wlg;
    
    // MPI
    mpi4(p,pgc,b->U);
    mpi4(p,pgc,UH);
    
    if(p->j_dir==1)
    {
    mpi4(p,pgc,b->V);
    mpi4(p,pgc,VH);
    }
    
    if(p->A220>0)
    {
    mpi4(p,pgc,b->W);
    mpi4(p,pgc,WH);
    }
    
    // physical boundaries: velocities
    GCSL4LOOP
    {
    i  = p->gcbsl4[n][0];
    j  = p->gcbsl4[n][1];
    cs = p->gcbsl4[n][3];
    bc = p->gcbsl4[n][4];
    
        // patch boundaries are set by patchBC
        if(bc>=100)
        continue;
        
        // inflow side
        if(cs==1 && (bc==1 || bc==6))
        {
            if(inflow_flag==0)
            {
                // relaxation zone: zero gradient, otherwise closed
                if(p->B98==2)
                {
                neumann(p,b->U,cs);
                neumann(p,b->V,cs);
                neumann(p,b->W,cs);
                }
                
                else
                {
                mirror(p,b->U,cs,-1.0);
                mirror(p,b->V,cs,1.0);
                mirror(p,b->W,cs,1.0);
                }
            }
        continue;
        }
        
        // outflow side
        if(cs==4 && (bc==2 || bc==7 || bc==8))
        {
            // active wave absorption (outflow_flag==1): ghost cells are set by ioflow
            
            if(outflow_flag==2)
            {
            neumann(p,b->U,cs);
            neumann(p,b->V,cs);
            neumann(p,b->W,cs);
            }
            
            if(outflow_flag==0)
            {
            mirror(p,b->U,cs,-1.0);
            mirror(p,b->V,cs,1.0);
            mirror(p,b->W,cs,1.0);
            }
        continue;
        }
        
        // walls, free slip
        if(cs==1 || cs==4)
        {
        mirror(p,b->U,cs,-1.0);
        mirror(p,b->V,cs,1.0);
        mirror(p,b->W,cs,1.0);
        }
        
        if(cs==2 || cs==3)
        {
        mirror(p,b->U,cs,1.0);
        mirror(p,b->V,cs,-1.0);
        mirror(p,b->W,cs,1.0);
        }
    }
    
    // physical boundaries: conserved variables from the ghost velocities and water depth
    GCSL4LOOP
    {
    i  = p->gcbsl4[n][0];
    j  = p->gcbsl4[n][1];
    cs = p->gcbsl4[n][3];
    
        for(q=1;q<=3;++q)
        {
        int ii=0,jj=0;
        
        if(cs==1)
        ii=-q;
        
        if(cs==4)
        ii=q;
        
        if(cs==2)
        jj=q;
        
        if(cs==3)
        jj=-q;
        
        wlg = MAX(b->eta(i+ii,j+jj) + b->depth(i+ii,j+jj), 0.0);
        
        WL(i+ii,j+jj) = wlg;
        UH(i+ii,j+jj) = b->U(i+ii,j+jj)*wlg;
        VH(i+ii,j+jj) = b->V(i+ii,j+jj)*wlg;
        WH(i+ii,j+jj) = b->W(i+ii,j+jj)*wlg;
        }
    }
}

void sflow_momentum_func::mirror(lexer *p, slice &f, int cs, double sgn)
{
    // symmetric (sgn=1) or antisymmetric (sgn=-1) extension across the boundary face
	if(cs==1)
	for(q=0;q<3;++q)
	f(i-q-1,j) = sgn*f(i+q,j);
	
	if(cs==2)
	for(q=0;q<3;++q)
	f(i,j+q+1) = sgn*f(i,j-q);
    
	if(cs==3)
	for(q=0;q<3;++q)
	f(i,j-q-1) = sgn*f(i,j+q);
	
	if(cs==4)
	for(q=0;q<3;++q)
	f(i+q+1,j) = sgn*f(i-q,j);
}

void sflow_momentum_func::neumann(lexer *p, slice &f, int cs)
{
	if(cs==1)
	for(q=0;q<3;++q)
	f(i-q-1,j) = f(i,j);
	
	if(cs==2)
	for(q=0;q<3;++q)
	f(i,j+q+1) = f(i,j);
    
	if(cs==3)
	for(q=0;q<3;++q)
	f(i,j-q-1) = f(i,j);
	
	if(cs==4)
	for(q=0;q<3;++q)
	f(i+q+1,j) = f(i,j);
}

void sflow_momentum_func::face_velocities(lexer *p, fdm2D *b, ghostcell *pgc, slice &WL)
{
    // staggered diagnostics for the sediment, turbulence, printing and ioflow modules
    SLICELOOP1
    {
    if(p->wet[IJ]==1 && p->wet[Ip1J]==1)
    b->P(i,j) = 0.5*(b->U(i,j) + b->U(i+1,j));
    
    else
    b->P(i,j) = 0.0;
    
    b->hx(i,j) = MAX(0.5*(WL(i,j) + WL(i+1,j)), 0.0);
    }
    
    SLICELOOP2
    {
    if(p->wet[IJ]==1 && p->wet[IJp1]==1)
    b->Q(i,j) = 0.5*(b->V(i,j) + b->V(i,j+1));
    
    else
    b->Q(i,j) = 0.0;
    
    b->hy(i,j) = MAX(0.5*(WL(i,j) + WL(i,j+1)), 0.0);
    }
    
    SLICELOOP4
    b->ws(i,j) = b->W(i,j);
    
    pgc->gcsl_start1(p,b->P,gcval_u);
    
    // in- and outflow faces (discharge diagnostics Qin2D, Qout2D)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    b->P(i-1,j) = b->U(i-1,j);
    }
    
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
    b->P(i,j) = (p->wet[IJ]==1)?0.5*(b->U(i,j) + b->U(i+1,j)):0.0;
    }

	pgc->gcsl_start2(p,b->Q,gcval_v);
    pgc->gcsl_start4(p,b->ws,gcval_w);
    pgc->gcsl_start1(p,b->hx,gcval_eta);
	pgc->gcsl_start2(p,b->hy,gcval_eta);
}

void sflow_momentum_func::irhs(lexer *p, fdm2D *b)
{
	SLICELOOP4
	{
	b->F(i,j) += b->Fext(i,j);
	b->Fext(i,j)=0.0;
	}
}

void sflow_momentum_func::jrhs(lexer *p, fdm2D *b)
{
    SLICELOOP4
	{
	b->G(i,j) += b->Gext(i,j);
	b->Gext(i,j)=0.0;
	}
}

void sflow_momentum_func::krhs(lexer *p, fdm2D *b)
{
    SLICELOOP4
	{
	b->H(i,j) += b->Hext(i,j);
	b->Hext(i,j)=0.0;
	}
}

void sflow_momentum_func::clearrhs(lexer *p, fdm2D *b)
{
	SLICELOOP4
	{
	b->F(i,j)=0.0;
    b->G(i,j)=0.0;
    b->H(i,j)=0.0;
	}
}
