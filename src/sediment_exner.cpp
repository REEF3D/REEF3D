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
#include"sediment_fdm.h"
#include"topo_relax.h"
#include"sediment_fou.h"
#include"sediment_cds.h"
#include"sediment_cds_hj.h"
#include"sediment_wenoflux.h"
#include"sediment_weno_hj.h"
#include"sflow_bicgstab.h"
#include"sediment_mixture.h"
#include<math.h>

sediment_exner::sediment_exner(lexer* p, ghostcell* pgc) : q0(p),xvec(p),rhsvec(p),M(p),qbx(p),qby(p),qbn(p),cxn(p),cyn(p)
{
    noneq_ini=0;

	if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;
    
    
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    
    Ls = p->S20;
    
    frac_k = -1;
    
    prelax = new topo_relax(p);
    
    if(p->S32==1)
    pdx = new sediment_fou(p);
    
    if(p->S32==2)
    pdx = new sediment_cds(p);
    
    if(p->S32==3)
    pdx = new sediment_cds_hj(p);
    
    if(p->S32==4)
    pdx = new sediment_wenoflux(p);
    
    if(p->S32==5)
    pdx = new sediment_weno_hj(p);
    
    psolv = new sflow_bicgstab(p,pgc);
}

sediment_exner::~sediment_exner()
{
}

void sediment_exner::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{   
    // multi-fraction bed
    if(s->pmix!=nullptr)
    {
    start_mixture(p,pgc,s);
    return;
    }
    
    // eq.
    if(p->S33==0)
    SEDSLICELOOP
    s->qb(i,j)=s->qbe(i,j);
    
    // non-eq.
    if(p->S33>0)
    non_equillibrium_solve(p,pgc,s); 
    
    qb_clear(p,s);
    
    pgc->gcsl_start4(p,s->qb,1);
    
    // suspended qs
    if(p->S62==2)
    susp_qs(p,pgc,s);
    
    // Exner
    if(p->S31==1)
    topovel1(p,pgc,s);
    
    if(p->S31==2)
    topovel2(p,pgc,s);
    
    if(p->S31==3)
    topovel3(p,pgc,s);
    
    if(p->S100>0)
	filter(p,pgc,s->vz,p->S100,p->S101);

	
    // Bedch
    timestep(p,pgc,s);
    

    SEDSLICELOOP
    s->dh(i,j) = p->dtsed*s->vz(i,j);
    
    // NHFLOW: deposit the suspended sediment of columns that fell dry
    dry_deposit(p,s);

	
	SEDSLICELOOP
    s->bedzh(i,j) += s->dh(i,j);

	pgc->gcsl_start4(p,s->bedzh,1);
}


void sediment_exner::start_RK(lexer* p, ghostcell* pgc, sediment_fdm *s)
{   
    // eq.
    if(p->S33==0)
    SEDSLICELOOP
    s->qb(i,j)=s->qbe(i,j);
    
    // non-eq.
    if(p->S33>0)
    non_equillibrium_solve(p,pgc,s); 
    
    qb_clear(p,s);
    
    pgc->gcsl_start4(p,s->qb,1);
    
    // suspended qs
    if(p->S62==2)
    susp_qs(p,pgc,s);
    
    // Exner
    if(p->S31==1)
    topovel1(p,pgc,s);
    
    if(p->S31==2)
    topovel2(p,pgc,s);
    
    if(p->S31==3)
    topovel3(p,pgc,s);
    
    if(p->S100>0)
	filter(p,pgc,s->vz,p->S100,p->S101);

	
    // Bedch
    timestep(p,pgc,s);

    //SEDSLICELOOP
    //s->dh(i,j) = p->dtsed*s->vz(i,j);
}
















void sediment_exner::start_mixture(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
    sediment_mixture *m = s->pmix;
    
    m->exner_begin(p,pgc,s);
    
    for(int q=0;q<m->nf;++q)
    {
        // qbe, shields parameters and grain size of fraction q
        m->load_fraction(p,pgc,s,q);
        
        SLICELOOP4
        s->qbe(i,j) = (*m->qbe_k[q])(i,j)*m->fac(i,j);
        
        pgc->gcsl_start4(p,s->qbe,1);
        
        // eq.
        if(p->S33==0)
        SEDSLICELOOP
        s->qb(i,j)=s->qbe(i,j);
        
        // non-eq., relaxation state per fraction
        if(p->S33>0)
        {
            SLICELOOP4
            qbn(i,j) = (*m->qbn_k[q])(i,j);
            
            noneq_ini = m->noneq_ini_k[q];
            
            non_equillibrium_solve(p,pgc,s);
            
            SLICELOOP4
            (*m->qbn_k[q])(i,j) = qbn(i,j);
            
            m->noneq_ini_k[q] = noneq_ini;
        }
        
        qb_clear(p,s);
        
        pgc->gcsl_start4(p,s->qb,1);
        
        // suspended load exchange distributed with the active layer composition
        frac_k = q;
        
        if(p->S62==2)
        susp_qs(p,pgc,s);
        
        SLICELOOP4
        s->vz(i,j) = 0.0;
        
        // Exner
        if(p->S31==1)
        topovel1(p,pgc,s);
        
        if(p->S31==2)
        topovel2(p,pgc,s);
        
        if(p->S31==3)
        topovel3(p,pgc,s);
        
        if(p->S100>0)
        filter(p,pgc,s->vz,p->S100,p->S101);
        
        SLICELOOP4
        {
        (*m->vz_k[q])(i,j) = s->vz(i,j);
        (*m->qb_k[q])(i,j) = s->qb(i,j);
        }
    }
    
    frac_k = -1;
    
    m->exner_end(p,pgc,s);
    
    // total
    SLICELOOP4
    {
    s->vz(i,j) = 0.0;
    s->qb(i,j) = 0.0;
    
        for(int q=0;q<m->nf;++q)
        {
        s->vz(i,j) += (*m->vz_k[q])(i,j);
        s->qb(i,j) += (*m->qb_k[q])(i,j);
        }
    }
    
    pgc->gcsl_start4(p,s->vz,1);
    pgc->gcsl_start4(p,s->qb,1);
    
    // Bedch
    timestep(p,pgc,s);
    
    SEDSLICELOOP
    s->dh(i,j) = p->dtsed*s->vz(i,j);
    
    // sorting: active layer and substrate
    for(int q=0;q<m->nf;++q)
    SLICELOOP4
    (*m->dh_k[q])(i,j) = p->dtsed*(*m->vz_k[q])(i,j);
    
    // NHFLOW: suspended sediment of columns that fell dry, distributed with the active layer
    // composition like the suspended exchange (susp_ED)
    if(p->A10==5 && p->S12>0)
    SEDSLICELOOP
    for(int q=0;q<m->nf;++q)
    (*m->dh_k[q])(i,j) += (*m->F[q])(i,j)*s->dryd(i,j)/(1.0-p->S24);
    
    dry_deposit(p,s);
    
    SEDSLICELOOP
    s->bedzh(i,j) += s->dh(i,j);
    
	pgc->gcsl_start4(p,s->bedzh,1);
    
    m->bedchange(p,pgc,s,m->dh_k);
}

void sediment_exner::qb_clear(lexer *p, sediment_fdm *s)
{
    // no bedload outside the sediment cells: s->qb is only written on SEDSLICELOOP, stale values
    // in cells that stopped being sediment cells were read as upwind neighbours (and, with a
    // multi-fraction bed, summed nf times per step)
    SLICEBASELOOP
    if(p->flagslice4[IJ]<0 || p->DFBED[IJ]<0)
    s->qb(i,j) = 0.0;
}

void sediment_exner::dry_deposit(lexer *p, sediment_fdm *s)
{
    // NHFLOW: the suspended sediment of columns that fell dry since the last bed update
    // (nhflow_suspended_IM1::drysave) is deposited as bed change dh. Columns without an
    // erodible bed (DFBED<0, solids) never exchange with the bed, their dryd is discarded
    // instead of growing without limit.
    if(p->A10!=5 || p->S12==0)
    return;
    
    SEDSLICELOOP
    {
    s->dh(i,j) += s->dryd(i,j)/(1.0-p->S24);
    s->dryd(i,j) = 0.0;
    }
    
    SLICEBASELOOP
    if(p->flagslice4[IJ]<0 || p->DFBED[IJ]<0)
    s->dryd(i,j) = 0.0;
}
