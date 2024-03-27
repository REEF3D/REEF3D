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

#include"6DOF_obj.h"
#include"lexer.h"
#include"momentum.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

void sixdof_obj::iniPosition_RBM(lexer *p, ghostcell *pgc)
{
    // Store initial position of triangles
    
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {        
            tri_x0[n][q] = tri_x[n][q] - c_(0);
            tri_y0[n][q] = tri_y[n][q] - c_(1);
            tri_z0[n][q] = tri_z[n][q] - c_(2);
        }
    }
	
	// Initial rotation

	if (p->X101==1)
	{	
        phi = p->X101_phi*(PI/180.0);
        theta = p->X101_theta*(PI/180.0);
        psi = p->X101_psi*(PI/180.0);	
	
		for (n=0; n<tricount; ++n)
		{
			rotation_tri
				(p,-phi,-theta,-psi,tri_x[n][0],tri_y[n][0],tri_z[n][0],c_(0),c_(1),c_(2));
			rotation_tri
				(p,-phi,-theta,-psi,tri_x[n][1],tri_y[n][1],tri_z[n][1],c_(0),c_(1),c_(2));
			rotation_tri
				(p,-phi,-theta,-psi,tri_x[n][2],tri_y[n][2],tri_z[n][2],c_(0),c_(1),c_(2));
		}

        // Rotate mooring end point
        if (p->X313==1) 
        {
            for (int line=0; line < p->mooring_count; line++)
            {
			    rotation_tri(p,-phi,-theta,-psi,p->X311_xe[line],p->X311_ye[line],p->X311_ze[line],c_(0),c_(1),c_(2));
            }
        }
	}
	

	// Initialise quaternions (Goldstein p. 604)
	e_(0) = 
		 cos(0.5*phi)*cos(0.5*theta)*cos(0.5*psi) 
		+ sin(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e_(1) = 
		 sin(0.5*phi)*cos(0.5*theta)*cos(0.5*psi) 
		- cos(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e_(2) = 
		 cos(0.5*phi)*sin(0.5*theta)*cos(0.5*psi) 
		+ sin(0.5*phi)*cos(0.5*theta)*sin(0.5*psi);
	e_(3) = 
		 cos(0.5*phi)*cos(0.5*theta)*sin(0.5*psi) 
		- sin(0.5*phi)*sin(0.5*theta)*cos(0.5*psi);   
        
        
    //cout<<"e_ "<<e_(0)<<" "<<e_(1)<<" "<<e_(2)<<" "<<e_(3)<<" "<<endl;
    
    en1_ = e_;
    en2_ = e_;
    en3_ = e_;
    ek_ = e_;
    
    cn1_ = c_;
    cn2_ = c_;
    cn3_ = c_;
    ck_ = c_;
    
    pn1_ = p_;
    pn2_ = p_;
    pn3_ = p_;
    pk_ = p_;
    
    hn1_ = h_;
    hn2_ = h_;
    hn3_ = h_;
    hk_ = h_;  


    // Initialise rotation matrices
    quat_matrices();
   
/* 
    cout<<p->mpirank<<" R_ "<<R_(0,0)<<" "<<R_(0,1)<<" "<<R_(0,2)<<" "<<R_(1,0)<<" "<<R_(1,2)<<" "<<R_(1,1)<<" "<<R_(2,0)<<" "<<R_(2,1)<<" "<<R_(2,2)<<" "<<endl;
    cout<<p->mpirank<<" G_ "<<G_(0,0)<<" "<<G_(0,1)<<" "<<G_(0,2)<<" "<<G_(1,0)<<" "<<G_(1,2)<<" "<<G_(1,1)<<" "<<G_(2,0)<<" "<<G_(2,1)<<" "<<G_(2,2)<<" "<<endl;
    cout<<p->mpirank<<" E_ "<<E_(0,0)<<" "<<E_(0,1)<<" "<<E_(0,2)<<" "<<E_(1,0)<<" "<<E_(1,2)<<" "<<E_(1,1)<<" "<<E_(2,0)<<" "<<E_(2,1)<<" "<<E_(2,2)<<" "<<endl;
    */
}

