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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"6DOF_obj_nhflow.h"
#include"6DOF_obj_cfd.h"
#include"6DOF_obj_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj_cfd::solve_eqmotion_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_cfd(p, a, pgc, alpha[iter], finalize);
    
    // load models (ship module) sample the CFD velocity
    sixdof_fluid_cfd fluid(p,a,pgc);
    pfluid = &fluid;
    
    update_forces(p);
    
    pfluid = nullptr;
    
    if(p->N40==2 || p->N40==12 || p->N40==22)
    rk2(p,pgc,iter);
    
    if(p->N40==3 || p->N40==13 || p->N40==23 || p->N40==33)
    rk3(p,pgc,iter);
   
    // low-storage RK3: FCLS3 (4, 24), RKLS3_df (14, also N40=13 with X10>0), RKLS3 (44)
    if(p->N40==4 || p->N40==14 || p->N40==24 || p->N40==44)
    rkls3(p,pgc,iter);
}

void sixdof_obj_nhflow::solve_eqmotion_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_nhflow(p, d, pgc, alpha[iter], finalize);

    // load models (ship module) sample the NHFLOW velocity
    sixdof_fluid_nhflow fluid(p,d,pgc);
    pfluid = &fluid;
    
    update_forces(p);
    
    pfluid = nullptr;
    
    // membranes (X 330) attached to the body: added-mass stabilisation of the partitioned coupling
    if(p->X330>0)
    membrane_stabilisation(p,iter);
    
    // porous floating body: linearly implicit drag
    if(p->X16==1)
    porous_damping_nhflow(p,iter);
    
    if(p->A510==2)
    rk2(p,pgc,iter);
    
    if(p->A510==3)
    rk3(p,pgc,iter);
}

void sixdof_obj_2D::solve_eqmotion_sflow(lexer *p, ghostcell *pgc, int iter, bool finalize)
{
    update_forces(p);
    
    if(p->A210==2)
    rk2(p,pgc,iter);
    
    if(p->A210==3)
    rk3(p,pgc,iter);
}

void sixdof_obj_nhflow::solve_eqmotion_oneway_nhflow(lexer *p, ghostcell *pgc, int iter, bool finalize)
{
    if(p->A510==2)
    rk2(p,pgc,iter);
    
    if(p->A510==3)
    rk3(p,pgc,iter);       
}

void sixdof_obj_2D::solve_eqmotion_oneway_sflow(lexer *p, ghostcell *pgc, int iter, bool finalize)
{
    if(p->A210==2)
    rk2(p,pgc,iter);
    
    if(p->A210==3)
    rk3(p,pgc,iter);       
}
    
void sixdof_obj::rk2(lexer *p, ghostcell *pgc, int iter)
{   
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rk2(iter,p->dt);
}

void sixdof_obj::rk3(lexer *p, ghostcell *pgc, int iter)
{   
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rk3(iter,p->dt);
}

void sixdof_obj::rkls3(lexer *p, ghostcell *pgc, int iter)
{
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rkls3(iter,gamma[iter],zeta[iter],p->dt);
}

void sixdof_obj::solve_eqmotion_oneway_onestep(lexer *p, ghostcell *pgc, bool finalize)
{
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.step_onestep(p->dt);
}

void sixdof_obj::rk4(lexer *p, ghostcell *pgc, int iter)
{
    // classical RK4, stage-synchronous with fnpf_RK4
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rk4(iter,p->dt);
}

