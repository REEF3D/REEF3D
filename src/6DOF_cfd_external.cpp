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

// CFD bodies and an external solver (6DOF_external.h; REEF3D::CORAL with fluid_body in coral.txt).
//
// CFD integrates the loads at the start of each stage (start_cfd, from the fields the stage starts
// with) and moves the body with them; the stages start and produce the bodies at (fractions of dt)
//   SSP-RK3 (N 40 3, 13, 23, 33):        start 0, 1, 1/2        bodies at 1, 1/2, 1
//   SSP-RK2 (N 40 2, 12, 22):            start 0, 1             bodies at 1, 1
//   low-storage RK3 (N 40 4, 14, 24, 44): start 0, 8/15, 2/3    bodies at 8/15, 2/3, 1
//
// Two-way (fluid_body, X 10 1), as FNPF: the loads of the first stage (t^n) are committed. With the
// predictor (default) the first stage takes a predicted pose (Taylor series from the external state at
// t^n), the loads of the second stage (at that pose) go over as provisional and the external solver
// steps to t^n + dt; without it, it steps at the first stage. The later stages take its motion within the
// step. CFD has no instantaneous added mass: light bodies need the artificial added mass of the external
// solver (fluid_body ... added_mass). The PTO (X 500) and the mooring (X 310) of two-way bodies are not
// used: they belong to the external model.
// One-way (fluid_body ... one_way, fluid_member): the kernel moves the body (fixed with X 11 0, or
// prescribed), its loads of the first stage go to the external solver, which then steps after the fluid
// step (one fluid step behind). Several bodies (X 20); not with the subcycling of the mesh refinement
// (G 7 1) yet.

#include"6DOF_cfd.h"
#include"6DOF_external.h"
#include"lexer.h"
#include"ghostcell.h"

bool sixdof_cfd::ext_drives(int nb) const
{
    return pext!=nullptr && pext->drives(nb) && !pext->one_way(nb);
}

void sixdof_cfd::ext_nodes(lexer *p, int iter, double &cs, double &cn, int &ns) const
{
    // start node of stage iter and time of the body it produces (fractions of dt), number of stages
    static const double c3[3] = {0.0, 1.0, 0.5}, n3[3] = {1.0, 0.5, 1.0};
    static const double c2[2] = {0.0, 1.0}, n2[2] = {1.0, 1.0};
    static const double cl[3] = {0.0, 8.0/15.0, 2.0/3.0}, nl[3] = {8.0/15.0, 2.0/3.0, 1.0};

    cs = 0.0;
    cn = 1.0;
    ns = 1;

    if(p->N40==3 || p->N40==13 || p->N40==23 || p->N40==33)
    {
        ns = 3;
        if(iter>=0 && iter<3) {cs = c3[iter]; cn = n3[iter];}
    }

    if(p->N40==2 || p->N40==12 || p->N40==22)
    {
        ns = 2;
        if(iter>=0 && iter<2) {cs = c2[iter]; cn = n2[iter];}
    }

    if(p->N40==4 || p->N40==14 || p->N40==24 || p->N40==44)
    {
        ns = 3;
        if(iter>=0 && iter<3) {cs = cl[iter]; cn = nl[iter];}
    }
}

void sixdof_cfd::ext_connect(lexer *p, ghostcell *pgc)
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

    if((twoway && p->X10!=1) || (p->G1>0 && p->G7==1))
    {
        if(p->mpirank==0)
        {
        if(twoway && p->X10!=1)
        cout<<"6DOF CFD: two-way external dynamics (fluid_body without one_way) needs a floating body (X 10 1)"<<endl;

        else
        cout<<"6DOF CFD: external dynamics (fluid_body, fluid_member) with the subcycling of the mesh refinement (G 7 1) is not supported yet"<<endl;
        }

        pgc->final();
        exit(1);
    }

    pext = e;
    ext_twoway = twoway;

    int n1 = 0;

    for(int nb=0; nb<number6DOF; ++nb)
    if(pext->drives(nb))
    {
        if(pext->one_way(nb))
        ++n1;

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

        if(p->X60!=1 && p->mpirank==0)
        cout<<"6DOF CFD: distributed loads (fluid_member) need the STL load integration (X 60 1)"<<endl;
        }
    }

    if(p->mpirank==0)
    {
        double cs, cn;
        int ns;
        ext_nodes(p,0,cs,cn,ns);

        if(ext_twoway)
        {
        cout<<"6DOF CFD: "<<number6DOF-n1<<(number6DOF-n1>1?" bodies follow":" body follows")<<" the external dynamics ("
        <<(pext->predictor() && ns>1 ? "predictor" : "no predictor")<<")"<<endl;

        if(p->X500>0)
        cout<<"6DOF CFD: the PTO X 500 is not used with external dynamics (PTO in the external model)"<<endl;

        if(p->X310>0)
        cout<<"6DOF CFD: the mooring X 310 is not used with external dynamics (mooring in the external model)"<<endl;
        }

        if(n1>0)
        cout<<"6DOF CFD: the loads of "<<n1<<(n1>1?" bodies go":" body goes")<<" one-way to the external solver"<<endl;
    }
}

void sixdof_cfd::ext_loads(lexer *p, double t, bool provisional)
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

void sixdof_cfd::ext_stage(lexer *p, int iter)
{
    // the loads of this stage were just integrated (start_cfd)
    const double t0 = p->simtime;
    const double t1 = p->simtime + p->dt;
    double cs, cn;
    int ns;
    ext_nodes(p,iter,cs,cn,ns);
    const bool pc = pext->predictor() && ns>1;

    // first stage: the committed loads of t^n (one-way bodies too)
    if(iter==0)
    {
        ext_loads(p,t0,false);

        if(!pc && ext_twoway)
        pext->advance(t0,t1);
    }

    // second stage with the predictor: its loads at the predicted pose, then the external step
    if(iter==1 && pc && ext_twoway)
    {
        ext_loads(p,t0 + cs*p->dt,true);
        pext->advance(t0,t1);
    }

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
