/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"momentum.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

void sixdof_df_object::ini_fbvel(lexer *p, fdm *a, ghostcell *pgc)
{
         double x0, x1, x2, y0, y1, y2, z0, z1, z2;
		double n0, n1, n2;
		double f1x,f2x,f3x,g0x,g1x,g2x,f1y,f2y,f3y;
		double g0y,g1y,g2y,f1z,f2z,f3z,g0z,g1z,g2z;
		double *integ;
		
		p->Darray(integ, 10);

		for (int n = 0; n < tricount; ++n)
		{
			x0 = tri_x[n][0];
			x1 = tri_x[n][1];
			x2 = tri_x[n][2];
		
			y0 = tri_y[n][0];
			y1 = tri_y[n][1];
			y2 = tri_y[n][2];
		
			z0 = tri_z[n][0];
			z1 = tri_z[n][1];
			z2 = tri_z[n][2];  
			
			n0 = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
			n1 = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
			n2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
		
			geometry_f(x0,x1,x2,f1x,f2x,f3x,g0x,g1x,g2x); 
			geometry_f(y0,y1,y2,f1y,f2y,f3y,g0y,g1y,g2y); 
			geometry_f(z0,z1,z2,f1z,f2z,f3z,g0z,g1z,g2z);	
		
			integ[0] += n0 * f1x;
			integ[1] += n0 * f2x;
			integ[2] += n1 * f2y;
			integ[3] += n2 * f2z;
			integ[4] += n0 * f3x;
			integ[5] += n1 * f3y;
			integ[6] += n2 * f3z;
			integ[7] += n0 * (y0 * g0x + y1 * g1x + y2 * g2x);
			integ[8] += n1 * (z0 * g0y + z1 * g1y + z2 * g2y);
			integ[9] += n2 * (x0 * g0z + x1 * g1z + x2 * g2z);	
		}

        double Vfb = integ[0]/6.0;
        double Rfb = 0.0;
        
		if (p->X22 == 1)
		{
			Mass_fb = p->X22_m;
			Rfb = Mass_fb/Vfb;
		}	
		else if (p->X21 == 1)
		{
			Rfb = p->X21_d;
			Mass_fb = Vfb*Rfb;
            p->X22_m = Mass_fb;
		}
        

    // Rigid body motion
    
    R_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    e_ << 0.0, 0.0, 0.0, 0.0;
    p_ << Uext*Mass_fb, Vext*Mass_fb, Wext*Mass_fb;
    c_ << 0.0, 0.0, 0.0;
    h_ << 0.0, 0.0, 0.0;
    
    omega_B << 0.0, 0.0, 0.0;
    omega_I << 0.0, 0.0, 0.0;
    
	if (p->X102 == 1)
	{
		p_(0) += p->X102_u*Mass_fb;
		p_(1) += p->X102_v*Mass_fb;
		p_(2) += p->X102_w*Mass_fb;
	} 
    
	if (p->X103 == 1)
	{
		h_(0) = p->X103_p;
		h_(1) = p->X103_q;
		h_(2) = p->X103_r;
	}  
	
    // Velocities
	p->ufb = p->vfb = p->wfb = 0.0;
	p->pfb = p->qfb = p->rfb = 0.0; 
	p->ufbi = p->vfbi = p->wfbi = 0.0;
	p->pfbi = p->qfbi = p->rfbi = 0.0; 
    
	if(p->X210==1)
	{
        p->ufbi = p->X210_u;
        p->vfbi = p->X210_v;
        p->wfbi = p->X210_w;
	}
	
	if(p->X211==1)
	{
        p->pfbi = p->X211_p;
        p->qfbi = p->X211_q;
        p->rfbi = p->X211_r;
	}

    p->ufbn = p->ufbi;
    p->vfbn = p->vfbi;
    p->wfbn = p->wfbi;
    p->pfbn = p->pfbi;   
    p->qfbn = p->qfbi;   
    p->rfbn = p->rfbi;
    
    p->del_Darray(integ, 10);	
}
  
void sixdof_df_object::ini_parameter(lexer *p, fdm *a, ghostcell *pgc)
{
        Uext = Vext = Wext = Pext = Qext = Rext = 0.0; 
    
    if (p->X210 == 1)
    {
        Uext = p->X210_u;
        Vext = p->X210_v;
        Wext = p->X210_w;
    }
    if (p->X211 == 1)
    {
        Pext = p->X211_p;
        Qext = p->X211_q;
        Rext = p->X211_r;
    }
    if (p->X221==1)
    {
        //motion_vec(p,a,pgc);
        cout<<"not implemented yet"<<endl;
    }
    
	
    // Velocities
	p->ufb = p->vfb = p->wfb = 0.0;
	p->pfb = p->qfb = p->rfb = 0.0; 
	p->ufbi = p->vfbi = p->wfbi = 0.0;
	p->pfbi = p->qfbi = p->rfbi = 0.0; 
    
	if(p->X210==1)
	{
        p->ufbi = p->X210_u;
        p->vfbi = p->X210_v;
        p->wfbi = p->X210_w;
	}
	
	if(p->X211==1)
	{
        p->pfbi = p->X211_p;
        p->qfbi = p->X211_q;
        p->rfbi = p->X211_r;
	}

    p->ufbn = p->ufbi;
    p->vfbn = p->vfbi;
    p->wfbn = p->wfbi;
    p->pfbn = p->pfbi;   
    p->qfbn = p->qfbi;   
    p->rfbn = p->rfbi;
    
    //p->del_Darray(integ, 10);	

    // Positions
    phi = theta = psi = 0.0;
    
    
    // Mass
    Mass_fb = 0.0;
    
    
    // Forces
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;
    Ffb_ << 0.0, 0.0, 0.0;
    Mfb_ << 0.0, 0.0, 0.0;
    
    // Printing
	printtime = 0.0;
    printtimenormal = 0.0;
    p->printcount_sixdof = 0;
}

