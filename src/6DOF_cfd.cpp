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

#include"6DOF_cfd.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ddweno_f_nug.h"
#include"field4a.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include<algorithm>

sixdof_cfd::sixdof_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = p->X20;
    
    // the fluid-momentum coupling (X 198) fits each body to the velocity in the region of the level set
    // of the body, which holds all bodies (a->fb): one body only so far
    if(p->X198>0 && number6DOF>1)
    {
        if(p->mpirank==0)
        cout<<"!!! X 198 (fluid-momentum coupling) works with one body (X 20 1) so far: X 197 / the explicit coupling instead"<<endl;
        
        p->X198 = 0;
    }
    
    // the one-step scheme with two projections (N 40 6): the bodies take the midpoint rule
    // (sixdof_rigidbody::stage_mid); the fluid-momentum coupling has no stages for it yet
    if(p->X198>0 && p->N40==6)
    {
        if(p->mpirank==0)
        cout<<"!!! X 198 (fluid-momentum coupling) has no stages for N 40 6 yet: X 197 / the explicit coupling instead"<<endl;
        
        p->X198 = 0;
    }
    
    // X 198 1 with the SSP Runge-Kutta stages (N 40 3 with the CFD mesh refinement, G 1 > 0; without it
    // the bodies take the low-storage RK3): the 2D swimmer is 9 % slow at N 47 0.25, as fast as with
    // the low-storage RK3 at 0.1 or with X 198 2 (refit after the projection; Fish Swimming validation 06)
    if(p->X198==1 && (p->N40==2 || p->N40==3 || p->N40==22 || p->N40==23 || p->N40==33) && p->N47>0.1 && p->mpirank==0)
    cout<<"!!! X 198 1 (fluid-momentum coupling) with the SSP Runge-Kutta stages: X 198 2 (the bodies refitted after each projection) or N 47 0.1 (at 0.25 a self-propelled swimmer is about 9 % slow)"<<endl;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj_cfd(p,pgc,nb));
    
    if(number6DOF>1)
    fball = new field4a(p);
    
    if(p->solidread>0 || p->topoforcing>0)
    {
    sf[0] = new field1(p);
    sf[1] = new field2(p);
    sf[2] = new field3(p);
    }
}
    
sixdof_cfd::~sixdof_cfd()
{
    delete fball;
    
    for(int q=0; q<3; ++q)
    delete sf[q];
}

void sixdof_cfd::start_cfd(lexer* p, fdm* a, ghostcell* pgc, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
    setup(p,a,pgc);
    
    // external dynamics (REEF3D::CORAL): the bodies are connected before their first stage
    if(!ext_checked)
    ext_connect(p,pgc);
    
    // the loads of all bodies first: the external solver takes them before it returns the poses
    if(pext!=nullptr && !amr_predict)
    {
        for (int nb=0; nb<number6DOF;++nb)
        fb_obj[nb]->hydrodynamic_forces_cfd(p,a,pgc,uvel,vvel,wvel,iter,finalize);
        
        ext_stage(p,iter);
    }
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Calculate forces (a predicted copy for the mesh refinement keeps the loads of the finest level)
        if(!amr_predict && pext==nullptr)
        fb_obj[nb]->hydrodynamic_forces_cfd(p,a,pgc,uvel,vvel,wvel,iter,finalize);
        
        // Advance body in time (two-way external bodies: the pose of the external solver, ext_stage)
        if(!ext_drives(nb))
        fb_obj[nb]->solve_eqmotion_cfd(p,a,pgc,iter,finalize);
         
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // deforming body (X 195): the deformation at the time of the stage
        fb_obj[nb]->deform_stage(p,iter);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_3D(p,a,pgc,finalize);  //----> main time consumer

        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        for(int q=0; q<3; ++q)
        fb_obj[nb]->sfh[q] = sf[q];
        
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        
        if(finalize)
        fb_obj[nb]->deform_step_end();
        
        for(int q=0; q<3; ++q)
        fb_obj[nb]->sfh[q] = nullptr;
        
        // several bodies: the level set of this one into that of all
        fb_combine(p,a,nb);
        
        // Print
        if(finalize==true && !amr_predict)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
            fb_obj[nb]->lateral_line(p,a,pgc);
        }
    }
    
    // several bodies: a->fb the level set of all bodies (ghost cell flags, implicit diffusion, level
    // set inside the bodies X 48 3)
    fb_finish(p,a,pgc);
    
    // actuator disks of the load models (ship propellers): momentum source of the fluid, after all
    // bodies have set their Heaviside fbh1..3 (no source in the forcing zone of any body)
    for (int nb=0; nb<number6DOF;++nb)
    fb_obj[nb]->actuator_forcing(p,a,pgc,fx,fy,fz);
    
    // ghostcell update
    pgc->gcdf_update(p,a);
}

void sixdof_cfd::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    
}

void sixdof_cfd::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}

// checkpoints (P 48, I 45): the bodies one after the other, each with its length
void sixdof_cfd::ckpt_state(lexer *p, std::vector<double> &v)
{
    v.push_back(double(number6DOF));
    for(int nb=0; nb<number6DOF; ++nb)
    {
        std::vector<double> s;
        fb_obj[nb]->ckpt_state(p,s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
}

bool sixdof_cfd::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
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
    return true;
}

// runtime coupling (reefcouple, I 243): CoG, velocity, angular velocity (inertial), rotation of body n
void sixdof_cfd::body_state(int n, double *s) const
{
    const sixdof_rigidbody &b = fb_obj[n]->rb;
    Eigen::Matrix<double,6,1> u;
    b.velocity(u);
    
    for(int q=0; q<3; ++q)
    {
    s[q]   = b.c(q);
    s[3+q] = u(q);
    s[6+q] = u(3+q);
    }
    
    for(int r=0; r<3; ++r)
    for(int q=0; q<3; ++q)
    s[9+3*r+q] = b.R(r,q);
}

// several bodies (X 20 > 1): the level set of all bodies, the minimum of their level sets
void sixdof_cfd::fb_combine(lexer *p, fdm *a, int nb)
{
    if(fball==nullptr)
    return;
    
    const int n = p->imax*p->jmax*p->kmax;
    
    if(nb==0)
    std::copy(a->fb.V, a->fb.V+n, fball->V);
    
    if(nb>0)
    for(int q=0; q<n; ++q)
    fball->V[q] = MIN(fball->V[q], a->fb.V[q]);
}

void sixdof_cfd::fb_finish(lexer *p, fdm *a, ghostcell *pgc)
{
    if(fball==nullptr)
    return;
    
    std::copy(fball->V, fball->V+p->imax*p->jmax*p->kmax, a->fb.V);
    
    pgc->start4a(p,a->fb,50);
}

// X 198 2: after the projection of stage iter, each body takes the rigid-body motion of the projected fluid in its
// hull region (the response of the added mass to this stage's forcing, which the fit of the next stage would see
// one stage later), and the fluid of the hull is forced again to it; the caller projects once more
void sixdof_cfd::refit_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, field &uvel, field &vvel, field &wvel,
                           field &fx, field &fy, field &fz, double alpha)
{
    for(int nb=0; nb<number6DOF; ++nb)
    fb_obj[nb]->momentum_refit(p,a,pgc,uvel,vvel,wvel);
    
    std::fill(fx.V, fx.V + p->imax*p->jmax*p->kmax, 0.0);
    std::fill(fy.V, fy.V + p->imax*p->jmax*p->kmax, 0.0);
    std::fill(fz.V, fz.V + p->imax*p->jmax*p->kmax, 0.0);
    
    // the weights of the hull of this stage again from those of the solids alone
    const int n = p->imax*p->jmax*p->kmax;
    if(sf[0]!=nullptr)
    {
    std::copy(sf[0]->V, sf[0]->V+n, a->fbh1.V);
    std::copy(sf[1]->V, sf[1]->V+n, a->fbh2.V);
    std::copy(sf[2]->V, sf[2]->V+n, a->fbh3.V);
    }
    else
    {
    std::fill(a->fbh1.V, a->fbh1.V+n, 0.0);
    std::fill(a->fbh2.V, a->fbh2.V+n, 0.0);
    std::fill(a->fbh3.V, a->fbh3.V+n, 0.0);
    }
    std::fill(a->fbh4.V, a->fbh4.V+n, 0.0);
    
    for(int nb=0; nb<number6DOF; ++nb)
    {
        for(int q=0; q<3; ++q)
        fb_obj[nb]->sfh[q] = sf[q];
        
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        
        for(int q=0; q<3; ++q)
        fb_obj[nb]->sfh[q] = nullptr;
    }
    
    ULOOP
    uvel(i,j,k) += alpha*p->dt*CPOR1*fx(i,j,k);
    VLOOP
    vvel(i,j,k) += alpha*p->dt*CPOR2*fy(i,j,k);
    WLOOP
    wvel(i,j,k) += alpha*p->dt*CPOR3*fz(i,j,k);
}

