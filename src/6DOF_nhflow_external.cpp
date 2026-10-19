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
Architect: Hans Bihs
--------------------------------------------------------------------*/

// NHFLOW bodies and an external solver (6DOF_external.h; REEF3D::CORAL with fluid_body in coral.txt).
//
// NHFLOW integrates the loads at the end of each stage (reforce_nhflow, after the projection); the 6DOF
// kernel moves the body at the start of the next stage with them. The stages produce the bodies at
//   RK3 (A 510 3, TVD): 1, 1/2, 1        RK2 (A 510 2): 1, 1        (fractions of dt)
//
// Two-way (fluid_body, X 10 1): instead of the kernel's RK stage the body takes the pose and velocity of
// the external solver. With its predictor (default) the first stage gets a predicted pose (Taylor series
// from the external state at t^n); the loads at the end of that stage go over as a provisional sample at
// t^n + dt and the external solver steps to t^n + dt, interpolating the loads. Without the predictor it
// steps at the start of the first stage, extrapolating from the committed samples. The later stages take
// the external solver's motion within the step. The loads at the end of the step are committed.
// NHFLOW has no instantaneous added mass: light bodies need the artificial added mass of the external
// solver (fluid_body ... added_mass).
// One-way (fluid_body ... one_way, fluid_member): the kernel moves the body (fixed with X 11 0, or
// prescribed), the loads at the end of the step go to the external solver.
// Both: the loads of the initialisation (X 10 1) are the first committed sample. The PTO (X 500) and
// the mooring (X 310) of two-way bodies are not used: they belong to the external model.
// X 13 (advance at the start of the stage, advance_stage): the same poses, set earlier.
// Mesh refinement with subcycling (G 7 1, nhflow_amr_sub.cpp): every step of the finest level is a fluid
// step of the external solver (stages and loads as above, amr_body_stage, amr_body_loads); the coarser
// levels, which step first, take poses predicted from the external state at the start of their step
// (amr_predict, the predicted copy of the 6DOF kernel). Without external bodies nothing changes.

#include"6DOF_nhflow.h"
#include"6DOF_external.h"
#include"lexer.h"
#include"ghostcell.h"

bool sixdof_nhflow::ext_drives(int nb) const
{
    return pext!=nullptr && pext->drives(nb) && !pext->one_way(nb);
}

double sixdof_nhflow::ext_node(lexer *p, int iter) const
{
    // time of the body state that stage iter produces (fraction of dt)
    static const double n3[3] = {1.0, 0.5, 1.0};

    if(p->A510==3 && iter>=0 && iter<3)
    return n3[iter];

    return 1.0;
}

void sixdof_nhflow::ext_connect(lexer *p, ghostcell *pgc)
{
    ext_checked = true;

    sixdof_external *e = sixdof_external::instance();

    if(e==nullptr)
    return;

    bool any = false, twoway = false;

    for(int nb=0; nb<number6DOF; ++nb)
    if(e->drives(nb))
    {
        any = true;

        if(!e->one_way(nb))
        twoway = true;
    }

    if(!any)
    return;

    if(twoway && p->X10!=1)
    {
        if(p->mpirank==0)
        cout<<"6DOF NHFLOW: two-way external dynamics (fluid_body without one_way) needs a floating body (X 10 1)"<<endl;

        pgc->final();
        exit(1);
    }

    pext = e;
    ext_twoway = twoway;

    for(int nb=0; nb<number6DOF; ++nb)
    if(pext->drives(nb))
    {
        double m;
        Eigen::Vector3d cg;
        Eigen::Matrix3d R, I;

        fb_obj[nb]->external_properties(m,cg,R,I);
        pext->connect(nb,m,cg,R,I);

        // distributed loads along a member of the external model
        Eigen::Vector3d a, ax;
        double L;
        int ns;

        if(pext->strips(nb,a,ax,L,ns))
        {
        fb_obj[nb]->external_strips_set(a,ax,L,ns);

        if(p->X60==2 && p->mpirank==0)
        cout<<"6DOF NHFLOW: distributed loads (fluid_member) need the STL load integration (X 60 1)"<<endl;
        }
    }

    if(p->mpirank==0)
    {
        if(ext_twoway)
        {
        cout<<"6DOF NHFLOW: the body follows the external dynamics ("<<(pext->predictor() ? "predictor" : "no predictor")<<")"<<endl;
        
        if(!pext->predictor())
        cout<<"6DOF NHFLOW: WARNING without the predictor the loads are extrapolated over the step; the NHFLOW loads react strongly to the body acceleration and the coupling can become unstable"<<endl;

        if(p->X500>0)
        cout<<"6DOF NHFLOW: the PTO X 500 is not used with external dynamics (PTO in the external model)"<<endl;

        if(p->X310>0)
        cout<<"6DOF NHFLOW: the mooring X 310 is not used with external dynamics (mooring in the external model)"<<endl;
        }

        else
        cout<<"6DOF NHFLOW: the loads of "<<number6DOF<<(number6DOF>1?" bodies go":" body goes")<<" one-way to the external solver"<<endl;
    }

    // the loads of the initialisation (floating bodies, X 10 1): the first committed sample (no strips yet);
    // not on a restart from a checkpoint (I 45 1), the external model has its samples
    if(p->X10!=1 || p->I45==1)
    return;
    
    std::vector<sixdof_external::vec6> F(number6DOF);
    std::vector<sixdof_external::mat6> A;

    for(int nb=0; nb<number6DOF; ++nb)
    fb_obj[nb]->external_fluid_load(p,F[nb]);

    pext->fluid_loads(p->simtime,F,A,false);
}

void sixdof_nhflow::ext_stage(lexer *p, int iter)
{
    const double t0 = p->simtime;
    const double t1 = p->simtime + p->dt;
    const bool pc = pext->predictor();

    if(iter==0 && !pc)
    pext->advance(t0,t1);

    const double cn = ext_node(p,iter);
    const double ts = (cn==1.0) ? t1 : t0 + cn*p->dt;

    for(int nb=0; nb<number6DOF; ++nb)
    if(ext_drives(nb))
    {
        Eigen::Vector3d xc, v, w;
        Eigen::Matrix3d R;

        if(iter==0 && pc)
        pext->predict(nb,ts,xc,R,v,w);

        else
        pext->state(nb,ts,xc,R,v,w);

        fb_obj[nb]->external_set_state(xc,R,v,w);
    }
}

void sixdof_nhflow::ext_loads(lexer *p, double t, bool provisional)
{
    std::vector<sixdof_external::vec6> F(number6DOF);
    std::vector<sixdof_external::mat6> A;

    for(int nb=0; nb<number6DOF; ++nb)
    {
    fb_obj[nb]->external_fluid_load(p,F[nb]);

    // the strips of the distributed loads (committed loads only)
    if(!provisional && fb_obj[nb]->external_strips_on())
    pext->fluid_strips(nb,t,fb_obj[nb]->external_strips());
    }

    pext->fluid_loads(t,F,A,provisional);
}

void sixdof_nhflow::ext_predict_copy(lexer *p, int iter)
{
    // a coarser level's copy of the bodies (subcycling): the pose predicted from the external state
    const double ts = p->simtime + ext_node(p,iter)*p->dt;
    
    for(int nb=0; nb<number6DOF; ++nb)
    if(ext_drives(nb))
    {
        Eigen::Vector3d xc, v, w;
        Eigen::Matrix3d R;
        
        pext->predict(nb,ts,xc,R,v,w);
        fb_obj[nb]->external_set_state(xc,R,v,w);
    }
}

void sixdof_nhflow::amr_body_stage(lexer *p, fdm_nhf *d, ghostcell *pgc, int s, bool predict)
{
    if(!ext_checked)
    ext_connect(p,pgc);
    
    if(ext_twoway)
    {
        if(predict)
        ext_predict_copy(p,s);
        
        else
        ext_stage(p,s);
    }
    
    for(int nb=0; nb<number6DOF; ++nb)
    fb_obj[nb]->amr_stage(p,d,pgc,s,predict,ext_drives(nb));
}

void sixdof_nhflow::amr_body_loads(lexer *p, ghostcell *pgc, int s, int ns)
{
    if(pext==nullptr)
    return;
    
    if(ext_twoway && s==0 && pext->predictor())
    {
        ext_loads(p,p->simtime + ext_node(p,0)*p->dt,true);
        pext->advance(p->simtime,p->simtime + p->dt);
    }
    
    if(s==ns-1)
    ext_loads(p,p->simtime + p->dt,false);
}
