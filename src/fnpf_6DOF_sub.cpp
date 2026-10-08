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


#include"fnpf_6DOF.h"
#include"6DOF_obj_fnpf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_fsf.h"
#include"fnpf_bed_update.h"
#include"fnpf_fsf_update.h"
#include"fnpf_amr.h"
#include"reefamr.h"
#include<functional>

//  Resolved bodies on a subcycled FNPF mesh refinement (fnpf_amr, G 7 1).
//
//  The finest level advances the body: in every stage of a finest step the loads from a psi_0
//  solve on the finest patches (the free-surface data from their tendencies, the parent columns
//  fixed at phi_t of the parent step), then the body stage with the added mass of the start of
//  the level-0 step.  Level 0 steps first: at its first stage the loads and the added mass from
//  the composite solves over all grids (all at the same time), then the stages of a predicted
//  copy (the loads frozen); the levels in between step with a predicted copy as well.  The
//  first stage of the first finest step takes the loads of level 0 (the same time).

void fnpf_6DOF::amr_save()
{
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->amr_save();
}

void fnpf_6DOF::amr_restore(lexer *p, ghostcell *pgc)
{
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->amr_restore(p,pgc);
}

// stage iter of the predicted copy with the step dt of its level
void fnpf_6DOF::amr_sub_predict(lexer *p, ghostcell *pgc, int iter, double dt)
{
    const double dts = p->dt;
    p->dt = dt;
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->amr_predict_fnpf(p,pgc,iter);
    
    p->dt = dts;
}

// stage iter of the finest level l (step dt from t): the loads (unless they are those of level 0),
// the body stage; the output at the last stage of the last finest step of the level-0 step
void fnpf_6DOF::amr_sub_finest(lexer *p, fdm_fnpf *c, ghostcell *pgc, int l, int iter, double t, double dt, bool loads, bool last)
{
    if(!initialized)
    return;
    
    const double dts = p->dt, ts = p->simtime;
    p->dt = dt;
    p->simtime = t;
    
    if(loads)
    forces_level(p,c,pgc,l);
    
    const bool finalize = (iter == ((p->A310==3) ? 2 : 3));
    
    for(int nb=0; nb<nbody; ++nb)
    {
        fb_obj[nb]->solve_eqmotion_fnpf(p,pgc,iter,finalize && last);
        fb_obj[nb]->update_position_fnpf(p,pgc,finalize && last);
        
        if(finalize && last)
        fb_obj[nb]->print_fnpf(p,pgc,iter);
    }
    
    p->dt = dts;
    p->simtime = ts;
}

// the body on the new sigma grids of the level-l patches (before their Laplace solve)
void fnpf_6DOF::amr_geometry_level(lexer *p, ghostcell *pgc, int l)
{
    if(!initialized)
    return;
    
    amr_grids(p,pgc);
    
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    if(amr->patch_level(G.id)==l)
    geometry(G,pgc);
}

// after it: wall nodes and body-band extrapolation of the level-l patches
void fnpf_6DOF::amr_post_solve_level(lexer *p, ghostcell *pgc, int l)
{
    if(!initialized)
    return;
    
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    if(amr->patch_level(G.id)==l)
    {
        amr->patch_walls_fi(G.id,G.c->Fi);
        extrapolate(G,pgc,G.c->Fi);
    }
}

// the loads of a stage of the finest level l: psi_0 on its patches (free-surface data from their
// tendencies, the parent columns at phi_t of the parent step), every hull triangle on the
// finest grid that holds its centroid; the added mass of the start of the level-0 step kept
void fnpf_6DOF::forces_level(lexer *p, fdm_fnpf *c, ghostcell *pgc, int l)
{
    amr_grids(p,pgc);
    
    const int ng = 1+(int)gp.size();
    vector<double*> fs(ng,nullptr);
    vector<slice*> Ds(ng,nullptr);
    fs[0] = g0.psi0;
    Ds[0] = &psiD;
    
    {
    reefamr_comms_off guard(pgc);
    for(size_t n=0; n<gp.size(); ++n)
    {
        fnpf_6DOF_grid &G = gp[n];
        fs[n+1] = G.psi0;
        Ds[n+1] = G.psiD;
        if(amr->patch_level(G.id)!=l)
        continue;
        
        // the velocities of the current state (the quadratic term of the pressure)
        G.pvel->velcalc_sig(G.p,G.c,pgc,G.c->Fi);
        body_velocities(G,pgc);
        
        lexer *p = G.p;
        slice4 &D = *G.psiD;
        slice &Ke = *G.Keta, &Kf = *G.Kfi;
        
        // psi_0: phi_t at fixed z on the free surface, phi_t|z = dFifsf/dt - Fz*deta/dt
        SLICELOOP4
        D(i,j) = Kf(i,j) - G.c->Fz(i,j)*Ke(i,j);
        
        pgc->gcsl_start4(p,D,50);
        
        zero_face(G);
        for(int nb=0; nb<nbody; ++nb)
        fb_obj[nb]->face_data_fnpf(G.p,G.c,pgc,-2,G.c->FBF,G.c->FBu,G.c->FBv,G.c->FBw);
        
        // Dirichlet on the free surface outside the footprint
        slice4 &foot = *G.foot;
        double *f = G.psi0;
        FFILOOP4
        if(foot(i,j)<0.5)
        {
            f[FIJK]   = D(i,j);
            f[FIJKp1] = D(i,j);
            f[FIJKp2] = D(i,j);
            f[FIJKp3] = D(i,j);
        }
        
        G.pbed->bedbc_sig(p,G.c,pgc,f,G.pf);
        amr->patch_walls_fi(G.id,f);
    }
    }
    
    // the columns around the patches: phi_t of the parent step
    amr->sub_psi_edges(l,&fs[0]);
    
    amr->lap_solve_psi_level(p,pgc,l,&fs[0],&Ds[0]);
    
    {
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    if(amr->patch_level(G.id)==l)
    {
        amr->patch_walls_fi(G.id,G.psi0);
        extrapolate(G,pgc,G.psi0);
    }
    }
    
    for(int nb=0; nb<nbody; ++nb)
    {
        sixdof_obj_fnpf::fnpf_force_sum S;
        fb_obj[nb]->forces_fnpf_zero(p,S);
        
        for(auto &G : gp)
        if(amr->patch_level(G.id)==l)
        {
            const int id = G.id;
            std::function<bool(double,double)> own = [&](double x, double y)
            {
                return amr->owns_point(x,y,id);
            };
            fb_obj[nb]->forces_fnpf_sum(G.p,G.c,G.psi0,G.psi,false,G.del,&own,S);
        }
        
        fb_obj[nb]->forces_fnpf_set(p,pgc,S,false);
    }
}
