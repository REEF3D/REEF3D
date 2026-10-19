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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_nhflow.h"
#include"6DOF_external.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
   
sixdof_nhflow::sixdof_nhflow(lexer *p, ghostcell *pgc) : press(p)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    // one body: d->FB, FHB and the pierced-column level (X 13) hold a single hull; every body's
    // ray cast would overwrite the previous one
    number6DOF = 1;
    
    if(p->X20>1)
    {
        if(p->mpirank==0)
        cout<<"6DOF NHFLOW: X 20 "<<p->X20<<" - only one floating body is supported in NHFLOW"<<endl;
        
        pgc->final();
        exit(1);
    }
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj_nhflow(p,pgc,nb));
}

sixdof_nhflow::~sixdof_nhflow()
{
}

void sixdof_nhflow::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, 
                                 double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    starttime = pgc->timer();
    
    // external dynamics (REEF3D::CORAL): the bodies are connected before their first stage
    if(!ext_checked)
    ext_connect(p,pgc);
    
    if(amr_predict && (p->X10==1 || p->X10==2))
    {
    if(ext_twoway && !advanced(p,iter))
    ext_predict_copy(p,iter);
    
    start_oneway(p,d,pgc,iter,FX,FY,FZ,WL,fe,false);
    p->fbtime+=pgc->timer()-starttime;
    return;
    }
    
    if(p->X10==1)
    start_twoway(p,d,pgc,iter,FX,FY,FZ,WL,fe,finalize);
    
    if(p->X10==2)
    start_oneway(p,d,pgc,iter,FX,FY,FZ,WL,fe,finalize);
    
    if(p->X10==3)
    start_shipwave(p,d,pgc,iter,finalize);
    
    if(p->X10==4)
    start_wavemaker(p,d,pgc,iter,FX,FY,FZ,WL,fe,finalize);
    
    p->fbtime+=pgc->timer()-starttime;
}

void sixdof_nhflow::start_twoway(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, 
                                double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    // two-way external dynamics: the poses of the stage come from the external solver
    if(ext_twoway && !advanced(p,iter))
    ext_stage(p,iter);
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Calculate forces
        //fb_obj[nb]->hydrodynamic_forces_nhflow(p,d,pgc,WL,finalize);
        
        // Advance body in time (X 13: already done at the start of the stage)
        if(!advanced(p,iter))
        {
        if(!ext_drives(nb))
        fb_obj[nb]->solve_eqmotion_nhflow(p,d,pgc,iter,finalize);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        }
        
        // Update position and trimesh
        fb_obj[nb]->update_position_nhflow(p,d,pgc,d->fs,finalize);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing_nhflow(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
        
        // Print
        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc,p->simtime+p->dt);
        }
    }
}

void sixdof_nhflow::start_oneway(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time (X 13: already done at the start of the stage)
        if(!advanced(p,iter))
        {
        if(!ext_drives(nb))
        fb_obj[nb]->solve_eqmotion_oneway_nhflow(p,pgc,iter,finalize);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        }
        
        // Update position and trimesh
        fb_obj[nb]->update_position_nhflow(p,d,pgc,d->fs,finalize);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing_nhflow(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
        
        // Print
        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            // Print
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc,p->simtime+p->dt);
        }
    }
}

void sixdof_nhflow::start_wavemaker(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    {

        // Update position and trimesh
        fb_obj[nb]->update_wavemaker_nhflow(p,d,pgc,d->fs,iter,finalize);  
        
        // Update forcing terms
        fb_obj[nb]->update_forcing_nhflow_wavemaker(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
        
        // Print
        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            // Print
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc,p->simtime+p->dt);
        }
    }
}

void sixdof_nhflow::start_shipwave(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, bool finalize)
{
    
    if(finalize==1)
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_oneway_onestep(p,pgc,finalize);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_2D(p,pgc,d->fs);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        if (p->X400==2)
        fb_obj[nb]->updateForcing_box(p,pgc,press);
        
        else if (p->X400==3)
        fb_obj[nb]->updateForcing_oned(p,pgc,press);
        
        else if (p->X400==10)
        fb_obj[nb]->updateForcing_stl(p,pgc,press,d->eta);
        
            // Print
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc,p->simtime+p->dt);
    }
}

void sixdof_nhflow::reforce_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, 
                                 double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    // shipwave (X 10 3): moving pressure field only, there is no body in the grid to re-impose
    // (d->FB is never cast, so the rigid-body forcing would act on the whole domain)
    if(p->X10==3)
    return;
    
    // subcycling (G 7 1): the predicted body of level 0 takes no loads
    if(amr_predict)
    {
    for (int nb=0; nb<number6DOF;++nb)
    {
        if(p->X10==4)
        fb_obj[nb]->update_forcing_nhflow_wavemaker(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
        
        else
        fb_obj[nb]->update_forcing_nhflow(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
    }
    return;
    }
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        // 1. re-impose the body velocity at the position the fluid was solved with:
        //    wavemaker (X 10 4): the paddle kinematics, not the rigid-body velocity u_fb
        //    porous body (X 16): the resistance is a momentum source and was applied once before
        //    the pressure projection - applying it again would double the drag. Only restore
        //    the porosity / FHB fields that reset() cleared.
        if(p->X10==4)
        fb_obj[nb]->update_forcing_nhflow_wavemaker(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
        
        else if(p->X16==0)
        fb_obj[nb]->update_forcing_nhflow(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
        
        else if(p->X16==1)
        fb_obj[nb]->porosity_nhflow(p,d,pgc);

        // 2. load from the pressure that was just solved
        fb_obj[nb]->hydrodynamic_forces_nhflow(p,d,pgc,WL,finalize);

        // 3. advance the body for the next substage
        /*fb_obj[nb]->solve_eqmotion_nhflow(p,d,pgc,iter,finalize);
        fb_obj[nb]->quat_matrices(p);
        fb_obj[nb]->update_position_nhflow(p,d,pgc,d->fs,finalize);
        fb_obj[nb]->update_fbvel(p,pgc);

        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            if(p->X50==1) fb_obj[nb]->print_vtp(p,pgc);
            if(p->X50==2) fb_obj[nb]->print_stl(p,pgc);
            fb_obj[nb]->print_parameter(p,pgc,p->simtime+p->dt);
        }*/
    }
    
    // external solver (CORAL): with the predictor the loads of the first stage at the predicted pose
    // (provisional), then the external step; the loads of the end of the step (committed)
    if(pext!=nullptr)
    {
        if(ext_twoway && iter==0 && pext->predictor())
        {
            ext_loads(p,p->simtime + ext_node(p,0)*p->dt,true);
            pext->advance(p->simtime,p->simtime + p->dt);
        }
        
        if(finalize)
        ext_loads(p,p->simtime + p->dt,false);
    }
}


void sixdof_nhflow::membrane_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, 
                                            double *UH, double *VH, double *WH, slice &WL)
{
    for (int nb=0; nb<number6DOF;++nb)
    fb_obj[nb]->membrane_forcing_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
}

void sixdof_nhflow::membrane_reaction_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, slice &WL, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    fb_obj[nb]->membrane_reaction_nhflow(p,d,pgc,alpha,WL,finalize);
}

bool sixdof_nhflow::membrane_iterated()
{
    bool it=false;
    
    for (int nb=0; nb<number6DOF;++nb)
    it = fb_obj[nb]->membrane_iterated() || it;
    
    return it;
}

void sixdof_nhflow::membrane_reforce_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    for (int nb=0; nb<number6DOF;++nb)
    fb_obj[nb]->membrane_reforce_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
}

bool sixdof_nhflow::membrane_couple_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double alpha, slice &WL, int it)
{
    // every body iterates its membranes; converged when all are
    bool conv=true;
    
    for (int nb=0; nb<number6DOF;++nb)
    conv = fb_obj[nb]->membrane_couple_nhflow(p,d,pgc,iter,alpha,WL,it) && conv;
    
    return conv;
}

// ------------------------------------------------------------------ X 13: advance before the level
bool sixdof_nhflow::advance_early_ok(lexer *p) const
{
    // the body ODE of a stage only needs the loads of the last stage. Nets (X 320), membranes
    // (X 330) and the ship load models (X 350) also sample the fluid of the stage in
    // solve_eqmotion_nhflow, so they keep the advance in phase_P1 (and the footprint of d->FB)
    return (p->X10==1 || p->X10==2) && p->X320==0 && p->X330==0 && p->X350==0 && p->X16==0;
}

bool sixdof_nhflow::advanced(lexer *p, int iter) const
{
    return adv_count==p->count && adv_iter==iter && adv_time==p->simtime;
}

void sixdof_nhflow::advance_stage(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter)
{
    if(advanced(p,iter))
    return;
    
    if(!ext_checked)
    ext_connect(p,pgc);
    
    const int last = (p->A510==2) ? 1 : 2;
    
    if(ext_twoway && p->X10==1)
    {
        if(amr_predict)
        ext_predict_copy(p,iter);
        
        else
        ext_stage(p,iter);
    }
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        if(ext_drives(nb))
        {}
        
        else if(amr_predict || p->X10==2)
        fb_obj[nb]->solve_eqmotion_oneway_nhflow(p,pgc,iter,false);
        
        else if(!ext_drives(nb))
        fb_obj[nb]->solve_eqmotion_nhflow(p,d,pgc,iter,iter==last);
        
        fb_obj[nb]->quat_matrices(p);
        fb_obj[nb]->pose_update(p,pgc);
    }
    
    adv_count = p->count;
    adv_iter = iter;
    adv_time = p->simtime;
}

bool sixdof_nhflow::top_phi_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, slice &WL, slice &phi)
{
    if(!advance_early_ok(p))
    {
        if(!adv_told && p->mpirank==0)
        cout<<"NHFLOW 6DOF X 13: footprint from the last ray cast (nets, membranes or ship loads sample the fluid in the body update)"<<endl;
        
        adv_told = true;
        return false;
    }
    
    advance_stage(p,d,pgc,iter);
    
    fb_obj[0]->top_phi_nhflow(p,d,WL,fb_obj[0]->nhflow_nb()*fb_obj[0]->nhflow_dsm(),phi);
    
    return true;
}

// checkpoints (P 48, I 45): the bodies one after the other, each with its length, and the stage the
// bodies were advanced for (X 13)
void sixdof_nhflow::ckpt_state(lexer *p, std::vector<double> &v)
{
    v.push_back(double(number6DOF));
    for(int nb=0; nb<number6DOF; ++nb)
    {
        std::vector<double> s;
        fb_obj[nb]->ckpt_state(p,s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
    v.push_back(double(adv_count)); v.push_back(double(adv_iter)); v.push_back(adv_time);
}

bool sixdof_nhflow::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
{
    size_t k=0;
    if(v.empty() || int(v[k++])!=number6DOF)
    return false;

    for(int nb=0; nb<number6DOF; ++nb)
    {
        const size_t n = (k<v.size()) ? size_t(v[k++]) : 0;
        if(n==0 || k+n>v.size())
        return false;
        if(!fb_obj[nb]->ckpt_restore(p,pgc,std::vector<double>(v.begin()+k,v.begin()+k+n)))
        return false;
        k += n;
    }
    if(k+3<=v.size())
    {
        adv_count = int(v[k++]); adv_iter = int(v[k++]); adv_time = v[k++];
    }
    return true;
}
