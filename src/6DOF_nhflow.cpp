/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_nhflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans.h"
   
sixdof_nhflow::sixdof_nhflow(lexer *p, ghostcell *pgc) : press(p)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = 1;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));
}

sixdof_nhflow::~sixdof_nhflow()
{
}

void sixdof_nhflow::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans* pvrans, vector<net*>& pnet, int iter, 
                                 double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    if(p->X10==1)
    start_twoway(p,d,pgc,pvrans,pnet,iter,FX,FY,FZ,WL,fe,finalize);
    
    if(p->X10==2)
    start_oneway(p,d,pgc,iter,FX,FY,FZ,WL,fe,finalize);
    
    if(p->X10==3)
    start_shipwave(p,d,pgc,finalize);
}

void sixdof_nhflow::start_twoway(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans* pvrans, vector<net*>& pnet, int iter, 
                                double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Calculate forces
        fb_obj[nb]->hydrodynamic_forces_nhflow(p,d,pgc,finalize);
        
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_nhflow(p,d,pgc,iter,pvrans,pnet);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices();
        
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
            
            fb_obj[nb]->print_parameter(p,pgc);
        }
    }
}

void sixdof_nhflow::start_oneway(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_oneway(p,pgc,iter);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices();
        
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
            
            fb_obj[nb]->print_parameter(p,pgc);
        }
    }
}

void sixdof_nhflow::start_shipwave(lexer *p, fdm_nhf *d, ghostcell *pgc, bool finalize)
{
    if(finalize==1)
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_oneway_onestep(p,pgc);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices();
        
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
        fb_obj[nb]->updateForcing_stl(p,pgc,press);
        
            // Print
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
    }
}
