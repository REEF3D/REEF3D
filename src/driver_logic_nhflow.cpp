/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include"vrans_header.h"
#include"nhflow_header.h"
#include"6DOF_void.h"
#include"6DOF_nhflow.h"

void driver::logic_nhflow()
{
    if(p->mpirank==0)
        cout<<"creating objects"<<endl;
    
    p->phimean = p->wd = p->F60;
    
    // nhflow
    if(p->A10!=5)
        pnhf = new nhflow_v(p,d,pgc);
    else if(p->A10==5)
        pnhf = new nhflow_f(p,d,pgc);
    
    // forcing
    pnhfdf = new nhflow_forcing(p);
    
    // FSF
    pnhfsf = new nhflow_fsf_f(p,d,pgc,pflow,pBC);
    
    // time stepping
    pnhfstep = new nhflow_timestep(p);

    //discretization scheme
    // signal speed
    pss = new nhflow_signal_speed(p);
    
    //Convection
    if(p->A511==1)
        pnhfconvec = new nhflow_HLL(p,pgc,pBC);
    else if(p->A511==2)
        pnhfconvec = new nhflow_HLLC(p,pgc,pBC);
    
    pnhfscalarconvec = new nhflow_scalar_ifou(p);
    
    //Diffusion
    if(p->A512==0)
        pnhfdiff = new nhflow_diff_void(p);
    else if(p->A512==1)
        pnhfdiff = new nhflow_ediff(p);
    else if(p->A512==2)
    {
        if( p->j_dir==0)
            pnhfdiff = new nhflow_idiff_2D(p);
        else if(p->j_dir==1)
            pnhfdiff = new nhflow_idiff(p);
    }
    
    // reconstruction
    if(p->A514<=3)
        precon = new nhflow_reconstruct_hires(p,pBC);
    else if(p->A514==4)
        precon = new nhflow_reconstruct_weno(p,pBC);
    else if(p->A514==5)
        precon = new nhflow_reconstruct_wenograd(p,pBC);
    
    //pressure scheme
    if(p->A520==0)
        pnhpress = new nhflow_pjm_hs(p,d,pBC);
    else if(p->A520==1)
        pnhpress = new nhflow_pjm(p,d,pgc,pBC);
    else if(p->A520==2)
        pnhpress = new nhflow_pjm_corr(p,d,pgc,pBC);

    //Turbulence
    if(p->A560==0)
        pnhfturb = new nhflow_komega_func_void(p,d,pgc);
    else if(p->A560==2 || p->A560==22)
    {
        pnhfturb = new nhflow_komega_IM1(p,d,pgc);
    
        if(p->j_dir==1)
            pnhfturbdiff = new nhflow_idiff(p);
        else if(p->j_dir==0)
            pnhfturbdiff = new nhflow_idiff_2D(p);
    }
    else if(p->A560==31)
        pnhfturb = new nhflow_LES_Smagorinsky(p,d,pgc);

    //Solver
    assign_solver();

    //Poison Solver
    assign_poisson_solver();
    
    //Data
    assign_data();
    
    //Printer
    if(p->P10==2)
        pnhfprint = new nhflow_vts3D(p,d,pgc);
    else
        pnhfprint = new nhflow_vtu3D(p,d,pgc);
    
    //VRANS
    assign_VRANS();
    
    //IOFlow
    assign_IOFlow();
    
    //Potential Flow Solver
    if(p->I11==0)
        pnhfpot = new nhflow_potential_v();
    else if(p->I11==1)
        pnhfpot = new nhflow_potential_f(p);
    
    //6DOF
    if(p->X10==0)
        p6dof = new sixdof_void(p,pgc);
    else if(p->X10>0)
        p6dof = new sixdof_nhflow(p,pgc);
    
    // Sediment
    if(p->S10==0)
        psed = new sediment_void();
    else if(p->S10>0)
        psed = new sediment_f(p,pgc,pturbcfd,pBC);
    
    //Momentum
    if(p->A510==2)
        pnhfmom = new nhflow_momentum_RK2(p,d,pgc,p6dof,pvrans,pnet,pnhfdf);
    else if(p->A510==3)
        pnhfmom = new nhflow_momentum_RK3(p,d,pgc,p6dof,pvrans,pnet,pnhfdf);
}
