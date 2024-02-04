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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_cfd.h"
#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ddweno_f_nug.h"

sixdof_cfd::sixdof_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = 1;
    
    for (int nb = 0; nb < number6DOF; nb++)
    p_df_obj.push_back(new sixdof_obj(p,a,pgc,nb));
}
    
sixdof_cfd::~sixdof_cfd()
{
}

	
void sixdof_cfd::start_forcing(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalise)
{
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Calculate forces
        p_df_obj[nb]->forces_stl(p,a,pgc,uvel,vvel,wvel,iter);
        
        // Advance body in time
        p_df_obj[nb]->solve_eqmotion(p,a,pgc,iter,pvrans,pnet);
        
        // Update position and fb level set
        p_df_obj[nb]->transform(p,a,pgc,finalise);  //----> main time consumer
        
        // Update forcing terms
        p_df_obj[nb]->updateForcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        
        // Save and print
        p_df_obj[nb]->interface(p,true);
        
        if(finalise==true)
        {
            p_df_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            p_df_obj[nb]->print_vtp(p,a,pgc);
            
            if(p->X50==2)
            p_df_obj[nb]->print_stl(p,a,pgc);
            
            p_df_obj[nb]->print_parameter(p, a, pgc);
        }
        
    }
    
    // ghostcell update
    pgc->gcdf_update(p,a);
}
