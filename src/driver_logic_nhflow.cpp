/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"waves_header.h"

void driver::logic_nhflow()
{    
	if(p->mpirank==0)
    cout<<"creating objects"<<endl;
    
    p->phimean = p->F60;
    
// nhflow
    if(p->A10!=55)
    pnh=new nhflow_v(p,d,pgc);
    
    if(p->A10==55)
    pnh=new nhflow_f(p,d,pgc);
    
    if(p->A10==55)
    {
    if(p->A540==1)
    pnhfsf = new nhflow_fsf_rk(p,d,pgc,pflow,pBC);
    
    if(p->A540==2)
    pnhfsf = new nhflow_fsf_fsm(p,d,pgc,pflow,pBC);
    }
    
// time stepping
    if(p->N48==0)
	ptstep=new fixtimestep(p);

    if((p->N48==1)  && (p->D20!=0&&p->D20!=2))
	ptstep=new etimestep(p);
	
	if((p->N48==1) && (p->D20==0||p->D20>=2))
	ptstep=new ietimestep(p);

//discretization scheme

    //Convection	
	/*if(p->D10==0)
	pconvec=new convection_void(p);

	if(p->D10==1)
	pconvec=new fou(p);

	if(p->D10==2)
	pconvec=new cds2(p);

	if(p->D10==3)
	pconvec=new quick(p);

	if(p->D10==4 && p->G2==0)
	pconvec=new weno_flux_nug(p);*/
    
    if(p->D10==4)
	pnhfconvec=new nhflow_weno_flux(p);
	
	/*if(p->D10==5)
	pconvec=new weno_hj_nug(p);
	
	if(p->D10==6)
	pconvec=new cds4(p);
    
    if(p->D10==7)
	pconvec=new weno3_flux(p);
    
    if(p->D10==8)
	pconvec=new weno3_hj(p);
	
	if(p->D10>=10 && p->D10<30)
	pconvec=new hires(p,p->D10);*/
    
    
//pressure scheme
	if(p->D30==0)
	ppress = new pressure_void(p);

    if(p->D30==1)
	ppress = new pjm_sig(p,a,pgc,pheat,pconc);
    
    if(p->D30==4)
	ppress = new pjm_sigss(p,a,pgc,pheat,pconc);

    if(p->D30==10)
	ppress = new pjm_sig_hs(p,a,pheat,pconc);


//poisson scheme for pressure
    if(p->D30==1)
	ppois = new poisson_sig(p,pheat,pconc);
	
//Solver
    if(p->j_dir==0)
	psolv = new bicgstab_ijk_2D(p,a,pgc);
    
    if(p->j_dir==1)
	psolv = new bicgstab_ijk(p,a,pgc);

//Poison Solver	
	if(p->N10==0)
	ppoissonsolv = new solver_void(p,a,pgc);
    
    if(p->N10==1 && p->j_dir==0)
	ppoissonsolv = new bicgstab_ijk_2D(p,a,pgc);
    
    if(p->N10==1 && p->j_dir==1)
	ppoissonsolv = new bicgstab_ijk(p,a,pgc);
	
	#ifdef HYPRE_COMPILATION
	if(p->N10>=10 && p->N10<20)
	ppoissonsolv = new hypre_struct(p,a,pgc,p->N10,p->N11);
	#endif
    
    #ifdef HYPRE_COMPILATION
	if(p->N10>=20 && p->N10<30)
	ppoissonsolv = new hypre_aij(p,a,pgc);
	#endif
    
    #ifdef HYPRE_COMPILATION
	if(p->N10>=30 && p->N10<40)
	ppoissonsolv = new hypre_sstruct(p,a,pgc);
	#endif
    
//IOFlow
	if(p->B60==0 && p->B90==0 && p->B180==0)
	pflow = new ioflow_v(p,pgc,pBC);

	if(p->B60>=1)
	pflow = new ioflow_f(p,pgc,pBC);

	if(p->B90>=1)
	pflow = new iowave(p,pgc,pBC);

	if(p->B180==1||p->B191==1||p->B192==1)
	pflow = new ioflow_gravity(p,pgc,pBC);
    
}