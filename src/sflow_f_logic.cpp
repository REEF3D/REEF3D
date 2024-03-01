/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"sflow_f.h"
#include"sflow_header.h"

void sflow_f::logic(lexer *p, fdm2D* b, ghostcell* pgc)
{	
    
	// timestep
    if(p->N48==0)
	ptime = new sflow_fixtimestep(p,b);
    
    if(p->N48==1)
	ptime = new sflow_etimestep(p,b);
	
	// convection
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
    
    if(p->A211==9)
    pconvec = new sflow_weno_blend(p);
    
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
        
        if(p->A211==9)
        pconvec = new sflow_weno_blend(p);
    }
    
    if(p->A215==1)
    {
        if(p->A211==0)
        pconvec = new sflow_voidconv(p);
        
        if(p->A211==1)
        pconvec = new sflow_cfou(p,b);
        
        if(p->A211==4 ||p->A211==5)
        pconvec = new sflow_cweno_flux(p,b);
        
        if(p->A211==6)
        pconvec = new sflow_chires(p,b,6);
        
        if(p->A211==7)
        pconvec = new sflow_chires(p,b,7);
        
        if(p->A211==8)
        pconvec = new sflow_chires(p,b,8);
        
        if(p->A211==9)
        pconvec = new sflow_weno_blend(p);
    }
 

    // filter
    pfilter = new sflow_filter(p);
	
	// free surface
	if(p->A240>=1)
	pfsf = new sflow_eta(p,b,pgc,pBC);

	// diffusion
	if(p->A212==0)
	pdiff =  new sflow_diffusion_void(p);
	
	if(p->A212==1)
	pdiff =  new sflow_ediff(p);
    
    if(p->A212==2)
	pdiff =  new sflow_idiff(p);
	
	// pressure
    if(p->A220==0)
	ppress = new sflow_hydrostatic(p,b,pBC);
    
    if(p->A220==1)
	ppress = new sflow_pjm_lin(p,b,pBC);
    
    if(p->A220==2)
	ppress = new sflow_pjm_quad(p,b,pBC);
    
    if(p->A220==3)
	ppress = new sflow_pjm_corr_lin(p,b,pBC);
    
    if(p->A220==5)
	ppress = new sflow_pjm_sw(p,b,pBC);
    
    // diffusion
	if(p->A260==0)
	pturb =  new sflow_turb_void(p);
    
    if(p->A260==1)
	pturb =  new sflow_turb_ke_IM1(p);
    
    if(p->A260==2)
	pturb =  new sflow_turb_kw_IM1(p);
    
    if(p->A260==3)
	pturb =  new sflow_turb_prandtl(p);
    
    if(p->A260==4)
	pturb =  new sflow_turb_parabolic(p);
    
    if(p->A260==5)
	pturb =  new sflow_turb_kw_IM1_v1(p);
    
    // Sediment
    if(p->S10==0)
    psed = new sediment_void();

    if(p->S10>0)
    psed = new sediment_f(p,aa,pgc,pturbcfd,pBC);
	
	// solver
	ppoissonsolv = new hypre_struct2D(p,pgc);
    
    
    psolv = new sflow_bicgstab(p,pgc);
    
    //IOFlow
	if(p->B60==0 && p->B90==0)
	pflow = new ioflow_v(p,pgc,pBC);

	if(p->B60>=1)
	pflow = new ioflow_f(p,pgc,pBC);

	if(p->B90>=1)
	pflow= new iowave(p,pgc,pBC);

	
	// printer
	pprint = new sflow_vtp_fsf(p,b,pgc);
	
	pprintbed = new sflow_vtp_bed(p,b);
    
    //6DOF
    if(p->X10!=3)
    p6dof = new sixdof_void(p,pgc);
    
    if(p->X10==3)
    p6dof = new sixdof_sflow(p,pgc);
	
	// momentum
    if(p->A210==1)
	pmom = new sflow_momentum_AB2(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,p6dof);
    
    if(p->A210==2)
	pmom = new sflow_momentum_RK2(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,p6dof);
    
	if(p->A210==3)
	pmom = new sflow_momentum_RK3(p,b,pconvec,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,p6dof);
    
    
    //Potential Flow Solver
    if(p->I11==0)
    potflow = new sflow_potential_v(p);

    if(p->I11==1)
    potflow = new sflow_potential_f(p);
    
}
