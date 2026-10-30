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
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"topo.h"
#include"reinitopo.h"
#include"suspended.h"
#include"bedload.h"
#include"bedshear.h"
#include"sandslide.h"
#include"topo_relax.h"
#include"bedslope.h"
#include"bedshear_reduction.h"
#include"bedload_direction.h"
#include"sediment_roughness.h"
#include"sediment_mixture.h"
#include"sediment_waves.h"
#include"bedload_waves.h"
#include"sediment_susp2D.h"
#include"slice1.h"
#include"slice2.h"

void sediment_f::sediment_algorithm_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &P, slice &Q)
{
    starttime=pgc->timer();
    
    ++p->sediter;
    
    // mesh refinement, multi-block grids: all grids together (sediment_amr.cpp)
    if(grids!=nullptr)
    {
    grids->sediment_step_sflow(p,b,pgc,pflow,this);
    }
    else
    {
    // prep, bed slope, roughness, bed shear stress (with the waves, S 110) -------
    physics_sflow(p,b,pgc,P,Q);

    // bedload *******
    if(s->pmix==nullptr)
    pbed->start(p,pgc,s);
    
    if(s->pmix!=nullptr)
    s->pmix->bedload_fractions(p,pgc,s,pbed);
    
    // wave formulas (S 11 6-9): the Exner step in the direction of the transport
    if(pbedw!=nullptr)
    transport_direction(p,pgc,1);
    
    // bedload_direction *******
    pbeddir->start(p,pgc,s);
	
    // relax *******
	prelax->start(p,pgc,s);
    
    // depth-averaged suspended load (S 120): the mean exchange rate since the last bed update
    if(psusp2D!=nullptr)
    psusp2D->bed_rate(p,pgc,s);
    
    // Exner *******
    ptopo->start(p,pgc,s);
    p->sedtime+=p->dtsed;
    
    if(psusp2D!=nullptr)
    psusp2D->bed_done(p);
    
    if(pbedw!=nullptr)
    transport_direction(p,pgc,0);
    
    // sandslide ********
    if(p->sediter%p->S94==0)
    pslide->start(p,pgc,s);
    
    // relax *******
	prelax->start(p,pgc,s);
    
    // update sflow  --------
    update_sflow(p,b,pgc,pflow);
    }
    
    // sediment print
    print_probes(p,pgc,s,pflow);
    
    // sediment log
    sedimentlog(p);
    
    if(s->pmix!=nullptr)
    s->pmix->print_log(p,pgc);
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Iter: "<<p->sediter<<" Sediment Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment CompTime: "<<setprecision(5)<<pgc->timer()-starttime<<endl;

}

// SFLOW: velocities, bed slope, roughness, bed shear stress; with the waves (S 110) the wave parameters and the
// wave and current bed shear stress before the critical bed shear stress
void sediment_f::physics_sflow(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    // prep SFLOW -------
    prep_sflow(p,b,pgc,P,Q);
    
    // the Eulerian current with the SEASTATE radiation stress (S 119 0)
    if(pwaves!=nullptr)
    pwaves->eulerian(p,b,pgc,s);
    
    // bedslope cds ******
    pslope->slope_cds(p,pgc,s);
    
    // bedslope reduction ******
    preduce->start(p,pgc,s);
    
    // bed roughness -------
    pks->start(p,pgc,s,b->hp);
    
    // bedshear stress -------
	pbedshear->taubed(p,b,pgc,s);
    
    if(pwaves!=nullptr)
    {
    pwaves->start(p,b,pgc,s);
    pwaves->bedshear(p,pgc,s);
    }
    
    pbedshear->taucritbed(p,b,pgc,s);
}

// SFLOW: depth-averaged suspended load (S 120), every flow step
void sediment_f::susp_sflow(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    physics_sflow(p,b,pgc,P,Q);
    
    // S 120 2: the suspended load of the transport formula
    if(p->S120==2 && pbedw!=nullptr)
    pbedw->start(p,pgc,s);
    
    psusp2D->start(p,b,pgc,s,pbedw);
}

// SFLOW wave formulas (S 11 6-9): the transport has its own direction (asymmetry, wave cycle); the Exner step,
// the bed slope correction of the direction (S 85) and the non-equilibrium relaxation (S 33) take the direction
// from s->P, s->Q, so these are replaced by the transport velocity |V| t during the step: t the unit transport
// direction (bedload_waves), |V| = max(|U|, Urms_eff), faces the mean of their two cells (closed next to solids),
// closed edge faces from the flow boundary conditions as in fill_PQ_sflow, open edge faces (a flow velocity there)
// from the zero-gradient ghost cells. Mode 0 restores the flow velocities.
void sediment_f::transport_direction(lexer *p, ghostcell *pgc, int mode)
{
    if(mode==0)
    {
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    s->P(i,j) = (*Psave)(i,j);
    s->Q(i,j) = (*Qsave)(i,j);
    }
    
    return;
    }
    
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    (*Psave)(i,j) = s->P(i,j);
    (*Qsave)(i,j) = s->Q(i,j);
    }
    
    SLICEBASELOOP
    (*Vtx)(i,j) = (*Vty)(i,j) = 0.0;
    
    SLICELOOP4
    {
    const double u = 0.5*(s->P(i,j)+s->P(i-1,j));
    const double v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
    double mag = sqrt(u*u + v*v);
    
    if(pwaves!=nullptr)
    mag = MAX(mag,pwaves->Ure(i,j));
    
    mag = MAX(mag,1.0e-6);
    
    (*Vtx)(i,j) = mag*pbedw->tx(i,j);
    (*Vty)(i,j) = mag*pbedw->ty(i,j)*p->y_dir;
    }
    
    pgc->gcsl_start4(p,*Vtx,1);
    pgc->gcsl_start4(p,*Vty,1);
    
    SLICELOOP1
    s->P(i,j) = (p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]>0) ? 0.5*((*Vtx)(i,j) + (*Vtx)(i+1,j)) : 0.0;
    
    SLICELOOP2
    s->Q(i,j) = (p->flagslice4[IJ]>0 && p->flagslice4[IJp1]>0) ? 0.5*((*Vty)(i,j) + (*Vty)(i,j+1)) : 0.0;
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);
    
    // open faces at the domain edges (in- and outflow: a flow velocity there): the transport velocity from the
    // zero-gradient ghost cells of Vtx, Vty; closed edges (walls) stay closed
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    if((i+p->origin_i==-1 || i+p->origin_i==p->gknox-1) && (*Psave)(i,j)!=0.0)
    s->P(i,j) = 0.5*((*Vtx)(i,j) + (*Vtx)(i+1,j));
    
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    if((j+p->origin_j==-1 || j+p->origin_j==p->gknoy-1) && (*Qsave)(i,j)!=0.0)
    s->Q(i,j) = 0.5*((*Vty)(i,j) + (*Vty)(i,j+1));
}
