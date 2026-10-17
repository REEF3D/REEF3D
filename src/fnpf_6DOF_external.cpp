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

// FNPF bodies driven by an external solver (6DOF_external.h; REEF3D::CORAL with fluid_body in
// coral.txt). The loads (psi_0, Bernoulli) and the added mass (unit-mode solves) are computed as for the
// 6DOF kernel; instead of the body RK stage they go to the external solver, which returns the pose and
// velocity at the time of the next stage:
//
//   RK3 (A 310 3, TVD): stage nodes 0, 1, 1/2, the stages produce the bodies at 1, 1/2, 1
//   RK4 (A 310 4):      stage nodes 0, 1/2, 1/2, 1, the stages produce the bodies at 1/2, 1/2, 1, 1
//
// First stage: the loads of the state at t^n (and the added mass) go to the external solver. Without
// its predictor it steps to t^n + dt at once. With the predictor (default) the next stage gets a
// predicted pose; its loads go over as provisional and the external solver steps then, interpolating
// the loads over the step. All later stages take the external solver's motion within the step.
// The PTO (X 500) and the mooring (X 310) of the 6DOF bodies are not used: they belong to the external
// model. X 11 selects the unit-mode solves of the added mass: the free DOFs there must contain every
// DOF the external joints let the body move in.

#include"fnpf_6DOF.h"
#include"6DOF_obj_fnpf.h"
#include"6DOF_external.h"
#include"fnpf_amr.h"
#include"lexer.h"
#include"ghostcell.h"

void fnpf_6DOF::ext_connect(lexer *p, ghostcell *pgc)
{
    ext_checked = true;
    
    sixdof_external *e = sixdof_external::instance();
    
    if(e==nullptr)
    return;
    
    int ndriven = 0;
    
    for(int nb=0; nb<nbody; ++nb)
    if(e->drives(nb))
    ++ndriven;
    
    if(ndriven==0)
    return;
    
    if(ndriven<nbody || (amr!=nullptr && (amr_on() || amr_sub())))
    {
        if(p->mpirank==0)
        {
        if(ndriven<nbody)
        cout<<"6DOF FNPF: "<<ndriven<<" of "<<nbody<<" bodies are driven externally (fluid_body), all or none must be"<<endl;
        else
        cout<<"6DOF FNPF: external dynamics (fluid_body) with mesh refinement (G 1) is not supported yet"<<endl;
        }
        
        pgc->final();
        exit(1);
    }
    
    pext = e;
    
    bool anyfixed = false;
    
    for(int nb=0; nb<nbody; ++nb)
    {
        double m;
        Eigen::Vector3d cg;
        Eigen::Matrix3d R, I;
        
        fb_obj[nb]->external_properties(m,cg,R,I);
        pext->connect(nb,m,cg,R,I);
        
        for(int n=0; n<6; ++n)
        if(fb_obj[nb]->dof_fixed_fnpf(p,n) && !(p->j_dir==0 && (n==1 || n==3 || n==5)))
        anyfixed = true;
    }
    
    if(p->mpirank==0)
    {
        cout<<"6DOF FNPF: "<<nbody<<(nbody>1?" bodies follow":" body follows")<<" the external dynamics ("
        <<(pext->predictor() ? "predictor" : "no predictor")<<")"<<endl;
        
        if(p->X500>0)
        cout<<"6DOF FNPF: the PTO X 500 is not used with external dynamics (PTO in the external model)"<<endl;
        
        if(p->X310>0)
        cout<<"6DOF FNPF: the mooring X 310 is not used with external dynamics (mooring in the external model)"<<endl;
        
        if(anyfixed)
        cout<<"6DOF FNPF: X 11 fixes DOFs: their added mass is not computed, the external joints must hold them"<<endl;
    }
}

void fnpf_6DOF::ext_stage_nodes(lexer *p, int iter, double &cs, double &cn) const
{
    // node of stage iter and time of the body state it produces (fractions of dt)
    static const double c3[3] = {0.0, 1.0, 0.5}, n3[3] = {1.0, 0.5, 1.0};
    static const double c4[4] = {0.0, 0.5, 0.5, 1.0}, n4[4] = {0.5, 0.5, 1.0, 1.0};
    
    cs = cn = 1.0;
    
    if(p->A310==3 && iter>=0 && iter<3)
    {cs = c3[iter]; cn = n3[iter];}
    
    if(p->A310!=3 && iter>=0 && iter<4)
    {cs = c4[iter]; cn = n4[iter];}
}

void fnpf_6DOF::ext_loads(lexer *p, double t, bool provisional, bool withA)
{
    std::vector<sixdof_external::vec6> F(nbody);
    std::vector<sixdof_external::mat6> A;
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->external_fluid_load(p,F[nb]);
    
    if(withA)
    {
        A.assign(nbody*nbody,sixdof_external::mat6::Zero());
        
        for(int nb=0; nb<nbody; ++nb)
        if(fb_obj[nb]->external_added_mass_on())
        A[nb*nbody+nb] = fb_obj[nb]->external_added_mass();
        
        // cross added mass, symmetrised as in motion_coupled
        if(nbody>1 && (int)Across.size()==nbody*nbody)
        for(int i=0; i<nbody; ++i)
        for(int j=0; j<nbody; ++j)
        if(i!=j)
        A[i*nbody+j] = 0.5*(Across[i*nbody+j] + Across[j*nbody+i].transpose());
        
        // DOFs without unit-mode solve (X 11): no added mass in their rows and columns
        for(int i=0; i<nbody; ++i)
        for(int n=0; n<6; ++n)
        if(fb_obj[i]->dof_fixed_fnpf(p,n))
        for(int j=0; j<nbody; ++j)
        {
        A[i*nbody+j].row(n).setZero();
        A[j*nbody+i].col(n).setZero();
        }
    }
    
    pext->fluid_loads(t,F,A,provisional);
}

void fnpf_6DOF::motion_external(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter)
{
    const int last = (p->A310==3) ? 2 : 3;
    const bool finalize = (iter==last);
    const double t0 = p->simtime;
    const double t1 = p->simtime + p->dt;
    const bool pc = pext->predictor();
    
    double cs, cn;
    ext_stage_nodes(p,iter,cs,cn);
    
    if(iter==0)
    {
        ext_loads(p,t0,false,true);
        
        if(!pc)
        pext->advance(t0,t1);
    }
    
    if(iter==1 && pc)
    {
        ext_loads(p,t0 + cs*p->dt,true,false);
        pext->advance(t0,t1);
    }
    
    const double ts = (cn==1.0) ? t1 : t0 + cn*p->dt;
    
    for(int nb=0; nb<nbody; ++nb)
    {
        Eigen::Vector3d xc, v, w;
        Eigen::Matrix3d R;
        
        if(iter==0 && pc)
        pext->predict(nb,ts,xc,R,v,w);
        
        else
        pext->state(nb,ts,xc,R,v,w);
        
        fb_obj[nb]->external_set_state(xc,R,v,w);
        fb_obj[nb]->update_position_fnpf(p,pgc,finalize);
        
        if(finalize)
        fb_obj[nb]->print_fnpf(p,pgc,iter);
    }
}
