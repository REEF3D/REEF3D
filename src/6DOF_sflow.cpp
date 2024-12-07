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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   
sixdof_sflow::sixdof_sflow(lexer *p, ghostcell *pgc) : press(p)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = 1;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));
}

sixdof_sflow::~sixdof_sflow()
{
}

void sixdof_sflow::start_sflow(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    if(p->X10==2)
    start_oneway(p,pgc,iter,fsglobal,P,Q,w,fx,fy,fz,finalize);
    
    if(p->X10==3)
    start_shipwave(p,pgc,iter,fsglobal,P,Q,fx,fy,fz,finalize);
}

void sixdof_sflow::start_oneway(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_oneway_sflow(p,pgc,iter);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_2D(p,pgc,fsglobal);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing_sflow(p,pgc,P,Q,w,fx,fy,fz,iter);
        
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

void sixdof_sflow::start_shipwave(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice&Q, slice &fx, slice &fy, slice &fz, bool finalize)
{
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        if(iter==0)
        fb_obj[nb]->solve_eqmotion_oneway_onestep(p,pgc);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_2D(p,pgc,fsglobal);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        if(p->X400==2)
        fb_obj[nb]->updateForcing_box(p,pgc,press);
        
        if(p->X400==3)
        fb_obj[nb]->updateForcing_oned(p,pgc,press);
        
        if(p->X400==10)
        fb_obj[nb]->updateForcing_stl(p,pgc,press);
        
            // Print
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
        
    }
}
