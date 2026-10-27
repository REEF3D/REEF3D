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

sixdof_cfd::sixdof_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = p->X20;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj_cfd(p,pgc,nb));
}
    
sixdof_cfd::~sixdof_cfd()
{
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
        
        // Update position and trimesh
        fb_obj[nb]->update_position_3D(p,a,pgc,finalize);  //----> main time consumer

        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        
        // Print
        if(finalize==true && !amr_predict)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
        }
    }
    
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
