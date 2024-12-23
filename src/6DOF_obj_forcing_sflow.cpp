/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::update_forcing_sflow(lexer *p, ghostcell *pgc, 
                             slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, int iter)
{
    // Calculate forcing fields
    double H, uf, vf, wf;
    double ef,efc;
    
    SLICELOOP4
    { 
        efc = 0.0;
        
        if(fs(i,j)<0.0)
        {
            efc = 0.0;
            
            if(fs(i-1,j)>0.0)   
            efc+=1.0;
            
            if(fs(i+1,j)>0.0)    
            efc+=1.0;

            if(fs(i,j-1)>0.0 && p->j_dir==1) 
            efc+=1.0;
            
            if(fs(i,j+1)>0.0 && p->j_dir==1)    
            efc+=1.0;
        }
        
        uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
        vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
        wf = 0.0;//u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
         
        if(efc>0.1)
        {
        H = Hsolidface_2D(p,1,0);
        fx(i,j) += H*(uf - P(i,j))/(alpha[iter]*p->dt);
        
        H = Hsolidface_2D(p,0,2);
        fy(i,j) += H*(vf - Q(i,j))/(alpha[iter]*p->dt);
        
        H = Hsolidface_2D(p,0,0);
        fz(i,j) += H*(wf - w(i,j))/(alpha[iter]*p->dt);
        }
    }
    
}
    
    
    
