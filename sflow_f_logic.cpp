/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"iowave.h"
#include"hypre_struct2D.h"
#include"sflow_etimestep.h"
#include"sflow_fou.h"
#include"sflow_cfou.h"
#include"sflow_weno_flux.h"
#include"sflow_weno_cflux.h"
#include"sflow_weno_hj.h"
#include"sflow_hires.h"
#include"sflow_chires.h"
#include"sflow_voidconv.h"
#include"sflow_eta.h"
#include"sflow_momentum_RK3.h"
#include"sflow_momentum_RK2.h"
#include"sflow_momentum_RK4.h"
#include"sflow_momentum_AB2.h"
#include"sflow_hydrostatic.h"
#include"sflow_vtp.h"
#include"sflow_vtp_bed.h"
#include"sflow_diffusion_void.h"
#include"sflow_ediff.h"
#include"sflow_pjm_lin.h"
#include"sflow_pjm_quad.h"
#include"sflow_pjm_sw.h"
#include"sflow_boussinesq_void.h"
#include"sflow_boussinesq_abbott.h"
#include"sflow_boussinesq_peregrine.h"
#include"sflow_boussinesq_madsen92.h"
#include"sflow_filter.h"

void sflow_f::logic(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	
	// timestep
	ptime = new sflow_etimestep(p,b);
	
	// convection
    if(p->A215==0)
    {
        if(p->A211==0)
        pconvec = new sflow_voidconv(p);
        
        if(p->A211==1)
        pconvec = new sflow_fou(p);
        
        if(p->A211==4)
        pconvec = new sflow_weno_flux(p);
        
        if(p->A211==5)
        pconvec = new sflow_weno_hj(p);
        
        if(p->A211==6)
        pconvec = new sflow_hires(p,6);
        
        if(p->A211==7)
        pconvec = new sflow_hires(p,7);
        
        if(p->A211==8)
        pconvec = new sflow_hires(p,8);
    }
    
    if(p->A215==1)
    {
        if(p->A211==0)
        pconvec = new sflow_voidconv(p);
        
        if(p->A211==1)
        pconvec = new sflow_cfou(p,b);
        
        if(p->A211==4 ||p->A211==5)
        pconvec = new sflow_weno_cflux(p,b);
        
        if(p->A211==6)
        pconvec = new sflow_chires(p,b,6);
        
        if(p->A211==7)
        pconvec = new sflow_chires(p,b,7);
        
        if(p->A211==8)
        pconvec = new sflow_chires(p,b,8);
    }
    
    // filter
    pfilter = new sflow_filter(p);
	
	// free surface
	if(p->A240>=1)
	pfsf = new sflow_eta(p,b,pgc);

	// diffusion
	if(p->A212==0)
	pdiff =  new sflow_diffusion_void(p);
	
	if(p->A212==1)
	pdiff =  new sflow_ediff(p);
	
	// pressure
    if(p->A220==0)
	ppress = new sflow_hydrostatic(p,b);
    
    if(p->A220==1)
	ppress = new sflow_pjm_lin(p,b);
    
    if(p->A220==2)
	ppress = new sflow_pjm_quad(p,b);
    
    if(p->A220==3)
	ppress = new sflow_pjm_sw(p,b);
    
    
    // Boussinesq wave model
    if(p->A230==0)
    pbouss = new sflow_boussinesq_void(p,b);
    
    if(p->A230==1)
    pbouss = new sflow_boussinesq_abbott(p,b);
    
    if(p->A230==2)
    pbouss = new sflow_boussinesq_peregrine(p,b);
    
    if(p->A230==3)
    pbouss = new sflow_boussinesq_madsen92(p,b);
	
	// solver
	ppoissonsolv = new hypre_struct2D(p,b,pgc);
    
	// ioflow
	pflow = new iowave(p,pgc);
	
	// printer
	pprint = new sflow_vtp(p,b,pgc);
	
	pprintbed = new sflow_vtp_bed(p,b);
	
	// momentum
	if(p->A210==3)
	pmom = new sflow_momentum_RK3(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,pbouss);
    
    if(p->A210==4)
	pmom = new sflow_momentum_RK4(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,pbouss);
	
    if(p->A210==11)
	pmom = new sflow_momentum_AB2(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf);
	
	if(p->A210==12)
	pmom = new sflow_momentum_RK2(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf);

	
    
    

}
