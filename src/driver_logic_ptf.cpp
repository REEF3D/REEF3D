/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"initialise_ptf.h"
#include"timestep_ptf.h"
#include"solver_ptf.h"
#include"bicgstab_ptf.h"
#include"bicgstab_ptf_2D.h"
#include"hypre_struct_ptf.h"
#include"solver_void_ptf.h"
#include"onephase_ptf.h"
#include"onephase_ptf_f.h"
#include"onephase_ptf_v.h"
void driver::logic_ptf()
{    
    if(p->mpirank==0)
    cout<<"creating objects"<<endl;
    
    pini = new initialise_ptf(p);

    if(p->mpirank==0)
	cout<<"starting ini"<<endl;
	pini->start(e,p,pgc);

// time stepping
    if(p->N48==0)
	ptstep=new fixtimestep_ptf(p);

	if(p->N48==1)
	ptstep=new pftimestep_ptf(p);
    
// Printer
	pprint = new vtu3D_ptf(p,e,pgc);
    
//IOFlow
	if(p->B60==0 && p->B90==0 && p->B180==0 )
	pflow = new ioflow_v(p,pgc,pBC);

	if(p->B90>=1)
	pflow= new iowave(p,pgc,pBC);
    
// Geodat
    if(p->G1==0)
    preto = new reinitopo_void();

    if(p->G1==1)
    {
    if(p->G40==0)
    preto = new reinitopo_void();
    
    if(p->G40==1)
    preto = new reinitopo_AB2(p);
    
    if(p->G40==3)
    preto = new reinitopo_RK3(p);
    }
    
//  Free Surface
    if(p->A10!=4)
    poneph_ptf = new onephase_ptf_v(p,e,pgc);
    
    if(p->A10==4)
    poneph_ptf = new onephase_ptf_f(p,e,pgc);
    
//  Laplace Solver	
	if(p->N10==0)
	plapsolv = new solver_void_ptf(p,e,pgc);
	
	if(p->N10==1)
	plapsolv = new bicgstab_ptf(p,e,pgc);
	
	#ifdef HYPRE_COMPILATION
	if(p->N10>10 && p->N10<=20)
    plapsolv = new hypre_struct_ptf(p,pgc,p->N10,p->N11);
    #endif
    
    #ifdef HYPRE_COMPILATION
	if(p->N10>20 && p->N10<=30)
	plapsolv = new hypre_aij_ptf(p,e,pgc);
	#endif
    
//  Voids
    
    pdata = new data_void_ptf(p,e,pgc);
    
    preini = new reini_void(p);
    
    pfsfdisc=new convection_void(p);
    
    pmp = new multiphase_v();
    
//  Wave Models
    if(p->A310==3)
    pptf = new ptf_RK3(p,e,pgc);
        
    if(p->A310==4)
    pptf = new ptf_RK4(p,e,pgc);
    
//Solver
    if(p->j_dir==0)
	psolv = new bicgstab_ptf_2D(p,e,pgc);

    if(p->j_dir==1)
	psolv = new bicgstab_ptf(p,e,pgc);
    
}
