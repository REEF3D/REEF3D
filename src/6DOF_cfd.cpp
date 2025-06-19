/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
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
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));
}
    
sixdof_cfd::~sixdof_cfd()
{
}

void sixdof_cfd::start_cfd(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
    setup(p,a,pgc);
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Calculate forces
        fb_obj[nb]->hydrodynamic_forces_cfd(p,a,pgc,uvel,vvel,wvel,iter,finalize);
        
        // Advance body in time
        fb_obj[nb]->solve_eqmotion(p,a,pgc,iter,pvrans,pnet);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_3D(p,a,pgc,finalize);  //----> main time consumer
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        
        
        // Print
        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
        }
    }
    
    // ghostcell update
    pgc->gcdf_update(p,a);
}

void sixdof_cfd::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    
}

void sixdof_cfd::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}
