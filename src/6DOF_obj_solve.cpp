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
#include"6DOF_motionext.h"

void sixdof_obj_cfd::solve_eqmotion_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_cfd(p, a, pgc, alpha[iter], finalize);
    
    // X 198: the fluid loads of the stage from the fitted velocity of the hull region
    if(p->X198>0)
    momentum_loads(p,iter);
    
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
    
    // the one-step scheme with two projections: position Verlet, stage 0 the half step (the loads
    // of time n), stage 1 the full step (the loads of the half step)
    if(p->N40==6)
    rkmid(p,pgc,iter);
    
    // the one-step scheme: one semi-implicit Euler step (the momenta with the loads of time n, the
    // pose with the new momenta)
    if(p->N40==5)
    {
    motion_stage_begin(p,1.0);
    
    get_trans(p,pgc);    
    get_rot(p);
    
    // the first step of a run from rest: there is no pressure yet (the loads are the weight alone),
    // the free DOFs would fall for a step and the fluid answer with a pressure pulse; they keep
    // their momenta, the prescribed DOFs (X 11 2) move
    if(p->count<=1 && p->I40==0)
    {
    rb.dp.setZero();
    rb.dh.setZero();
    }
    
    rb.step_symplectic(p->dt);
    
    motion_stage_end(p);
    }
}

void sixdof_obj_nhflow::solve_eqmotion_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_nhflow(p, d, pgc, alpha[iter], finalize);
    
    // power take-off (X 500): load of the stage, as in FNPF
    if(pto_on_)
    pto_forces(p,pgc,iter);

    // load models (ship module) sample the NHFLOW velocity
    sixdof_fluid_nhflow fluid(p,d,pgc);
    pfluid = &fluid;
    
    update_forces(p);
    
    pfluid = nullptr;
    
    // membranes (X 330) attached to the body: added-mass stabilisation of the partitioned coupling
    if(p->X330>0)
    membrane_stabilisation(p,iter);
    
    // the hull as a sharp thin body (X 15 8) or a solid with cut faces (X 15 10): added-mass stabilisation of the
    // partitioned coupling
    if((p->X15==8 || p->X15==10) && d->hullbody!=nullptr)
    hull_stabilisation(p,iter);
    
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
    // SSP-RK2: both stages end at t+dt
    motion_stage_begin(p,1.0);
    
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rk2(iter,p->dt);
    
    motion_stage_end(p);
}

void sixdof_obj::rk3(lexer *p, ghostcell *pgc, int iter)
{   
    // SSP-RK3 (Shu-Osher): the stages end at t+dt, t+dt/2, t+dt
    motion_stage_begin(p, iter==1 ? 0.5 : 1.0);
    
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rk3(iter,p->dt);
    
    motion_stage_end(p);
}

// The body with the one-step scheme with two projections (CFD N 40 6): position Verlet
// (Stoermer-Verlet, second order, symplectic: a floating body in heave keeps its energy).
//   iter 0, the half step: the pose from n to n+1/2 with the momenta of n, the momenta to n+1/2
//           with the loads of p(n) (the velocity the half step forces)
//   iter 1, the full step: the momenta from n to n+1 with the loads of p(n+1/2), the pose from
//           n+1/2 to n+1 with the momenta of n+1
// Prescribed motions (sixdof_motionext) act on the derivatives as in the stages: get_trans and
// get_rot again after the momenta, so the pose takes the prescribed velocity.
void sixdof_obj::rkmid(lexer *p, ghostcell *pgc, int iter)
{
    const double dt = p->dt;
    
    // the half step ends at t+dt/2, the full step at t+dt
    motion_stage_begin(p, iter==0 ? 0.5 : 1.0);
    
    get_trans(p,pgc);    
    get_rot(p);
    
    if(iter==0)
    {
        rb.pk = rb.p;
        rb.ck = rb.c;
        rb.hk = rb.h;
        rb.ek = rb.e;
        
        rb.c = rb.ck + 0.5*dt*rb.dc;
        rb.e = (rb.ek + 0.5*dt*rb.de).normalized();
        
        rb.p = rb.pk + 0.5*dt*rb.dp;
        rb.h = rb.hk + 0.5*dt*rb.dh;
        
        motion_stage_end(p);
        return;
    }
    
    rb.p = rb.pk + dt*rb.dp;
    rb.h = rb.hk + dt*rb.dh;
    
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.c = rb.c + 0.5*dt*rb.dc;
    rb.e = (rb.e + 0.5*dt*rb.de).normalized();
    
    motion_stage_end(p);
}

void sixdof_obj::rkls3(lexer *p, ghostcell *pgc, int iter)
{
    // low-storage RK3 (Williamson): the stages end at t+8/15 dt, t+2/3 dt, t+dt
    motion_stage_begin(p, iter==0 ? 8.0/15.0 : (iter==1 ? 2.0/3.0 : 1.0));
    
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rkls3(iter,gamma[iter],zeta[iter],p->dt);
    
    motion_stage_end(p);
}

void sixdof_obj::solve_eqmotion_oneway_onestep(lexer *p, ghostcell *pgc, bool finalize)
{
    // one step to t+dt (the SFLOW semi-implicit step sflow_momentum_SI, the ship waves X 10 3)
    motion_stage_begin(p,1.0);
    
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.step_onestep(p->dt);
    
    motion_stage_end(p);
}

void sixdof_obj::rk4(lexer *p, ghostcell *pgc, int iter)
{
    // classical RK4, stage-synchronous with fnpf_RK4: the stages end at t+dt/2, t+dt/2, t+dt, t+dt
    motion_stage_begin(p, iter<2 ? 0.5 : 1.0);
    
    get_trans(p,pgc);    
    get_rot(p);
    
    rb.stage_rk4(iter,p->dt);
    
    motion_stage_end(p);
}

// Prescribed motion from a file (X 240 1): the velocity and the position of the body at the time the stage of
// its integrator ends, the time of the fluid state that sees the body next.  The body moves before the fluid
// stage (CFD, NHFLOW, SFLOW: the forcing of the stage; FNPF: the body boundary of the next Laplace solve).
// The velocity of the file interval at p->simtime, the start of the step, in every stage was a step function
// in time and the same target for stages that end at different times.  Other motions ignore the time.
void sixdof_obj::motion_stage_begin(lexer *p, double cend)
{
    if(pmotion!=nullptr)
    pmotion->stage_begin(p, p->simtime + cend*p->dt, rb.c);
}

void sixdof_obj::motion_stage_end(lexer *p)
{
    if(pmotion!=nullptr)
    pmotion->stage_position(p,rb.c);
}

