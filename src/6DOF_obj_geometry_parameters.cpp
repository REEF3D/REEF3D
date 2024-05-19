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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::geometry_parameters(lexer *p, fdm *a, ghostcell *pgc)
{
    double x0, x1, x2, y0, y1, y2, z0, z1, z2;
    double n0, n1, n2;
    double f1x,f2x,f3x,g0x,g1x,g2x,f1y,f2y,f3y;
    double g0y,g1y,g2y,f1z,f2z,f3z,g0z,g1z,g2z;
    double *integ;
    
    I_ = Eigen::Matrix3d::Zero();

	if (p->X131 > 0 || p->X132 > 0 || p->X133 > 0 || p->X153 > 0)
	{
		geometry_ls(p,a,pgc);
	}	
        
	else
	{

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
        
        double Vol,Vol_ls,H;
        Vol=0.0;
        for (int n = 0; n < tricount; ++n)
        {
        double x1, x2, x3, y1, y2, y3, z1, z2, z3;
        
             x1 = tri_x[n][0];
			x2 = tri_x[n][1];
			x3 = tri_x[n][2];
		
			y1 = tri_y[n][0];
			y2 = tri_y[n][1];
			y3 = tri_y[n][2];
		
			z1 = tri_z[n][0];
			z2 = tri_z[n][1];
			z3 = tri_z[n][2];  
            
            Vol += (1.0/6.0)*(-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
            
        }
        
        Vol_ls=0.0;
        ALOOP
        {
            H = Hsolidface(p,a,0,0,0);
            Vol_ls+= p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H;
        }
        Vol_ls=pgc->globalsum(Vol_ls);

        if(p->X180==0)
        {
            Vfb = integ[0]/6.0;
            Rfb = 0.0;
            
            if (p->X22==1)
            {
                Mass_fb = p->X22_m;
                Rfb = Mass_fb/Vfb;
            }	
            
            else if (p->X21==1)
            {
                Rfb = p->X21_d;
                Mass_fb = Vfb*Rfb;
                p->X22_m = Mass_fb;
            }
        }

		integ[1] *= Rfb/24.0;
		integ[2] *= Rfb/24.0;
		integ[3] *= Rfb/24.0;

		c_(0) = integ[1]/Mass_fb;
		c_(1) = integ[2]/Mass_fb;
		c_(2) = integ[3]/Mass_fb;
		
		if(p->X23==1)
		{
			c_(0) = p->X23_x; 
			c_(1) = p->X23_y; 
			c_(2) = p->X23_z; 
		}
		
        double Ix, Iy, Iz;

		if(p->X24==1)
		{
			Ix = p->X24_Ix; 
			Iy = p->X24_Iy; 
			Iz = p->X24_Iz; 

			I_(0,1) = 0.0;
			I_(0,2) = 0.0;
			I_(1,2) = 0.0;      		
		}
		else
		{
			integ[4] *= Rfb/60.0;
			integ[5] *= Rfb/60.0;
			integ[6] *= Rfb/60.0;
			integ[7] *= Rfb/120.0;
			integ[8] *= Rfb/120.0;
			integ[9] *= Rfb/120.0;
		
			Ix = integ[5] + integ[6] - Mass_fb*(c_(1)*c_(1) + c_(2)*c_(2));
			Iy = integ[4] + integ[6] - Mass_fb*(c_(2)*c_(2) + c_(0)*c_(0));
			Iz = integ[4] + integ[5] - Mass_fb*(c_(0)*c_(0) + c_(1)*c_(1));
			I_(0,1) = -integ[7] + Mass_fb*c_(0)*c_(1);
			I_(1,2) = -integ[8] + Mass_fb*c_(1)*c_(2);
			I_(0,2) = -integ[9] + Mass_fb*c_(2)*c_(0);
		}

		I_(0,0) = Ix;
		I_(1,0) = I_(0,1);
		I_(1,1) = Iy;
		I_(2,0) = I_(0,2); 
		I_(2,1) = I_(1,2); 
		I_(2,2) = Iz;  
        p->W_fb = Rfb;

        p->xg = c_(0);
        p->yg = c_(1);
        p->zg = c_(2);

        p_ *= Mass_fb;

		if(p->mpirank==0)
		{
			cout<<"Center of Gravity xg: "<<c_(0)<<" yg: "<<c_(1)<<" zg: "<<c_(2)<<endl;
			cout<<"Volume Floating Body: "<<Vfb<<" Vol_ls: "<<Vol_ls<<endl;
			cout<<"Mass Floating Body: "<<Mass_fb<<endl;
			cout<<"Density Floating Body: "<<Rfb<<endl;
			cout<<"Moments of Inertia Tensor:\n"<<I_<<endl;
		}
		
		p->del_Darray(integ, 10);	
	}
}


void sixdof_obj::geometry_parameters_2D(lexer *p, ghostcell *pgc)
{
    double x0, x1, x2, y0, y1, y2, z0, z1, z2;
    double xmin = .10e9;
    double xmax =-.10e9;
    double ymin = .10e9;
    double ymax =-.10e9;
    double zmin = .10e9;
    double zmax =-.10e9;
    
    Mass_fb = 1.0;
    Vfb = 1.0;
    Rfb = Mass_fb/Vfb;
    
    I_ = Eigen::Matrix3d::Zero();

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
            
            
            xmin = MIN3(x0,x1,x2);
            xmax = MAX3(x0,x1,x2);
            
            ymin = MIN3(y0,y1,y2);
            ymax = MAX3(z0,z1,z2);
            
            zmin = MIN3(z0,z1,z2);
            zmax = MAX3(z0,z1,z2);
		}
        
        if(p->X23!=1)
		{
			c_(0) = xmin + 0.5*(xmax-xmin); 
			c_(1) = ymin + 0.5*(ymax-ymin); 
			c_(2) = zmin + 0.5*(zmax-zmin); 
		}

		if(p->X23==1)
		{
			c_(0) = p->X23_x; 
			c_(1) = p->X23_y; 
			c_(2) = p->X23_z; 
		}
		
        double Ix, Iy, Iz;

		if(p->X24!=1)
		{
			Ix = 1.0; 
			Iy = 1.0; 
			Iz = 1.0; 

			I_(0,1) = 0.0;
			I_(0,2) = 0.0;
			I_(1,2) = 0.0;      		
		}
        
        if(p->X24==1)
		{
			Ix = p->X24_Ix; 
			Iy = p->X24_Iy; 
			Iz = p->X24_Iz; 

			I_(0,1) = 0.0;
			I_(0,2) = 0.0;
			I_(1,2) = 0.0;      		
		}


		I_(0,0) = Ix;
		I_(1,0) = I_(0,1);
		I_(1,1) = Iy;
		I_(2,0) = I_(0,2); 
		I_(2,1) = I_(1,2); 
		I_(2,2) = Iz; 
 

        p->xg = c_(0);
        p->yg = c_(1);
        p->zg = c_(2);

        p_ *= Mass_fb;

		if(p->mpirank==0)
		{
			cout<<"Center of Gravity xg: "<<c_(0)<<" yg: "<<c_(1)<<" zg: "<<c_(2)<<endl;
			cout<<"Volume Floating Body: "<<Vfb<<endl;
			cout<<"Mass Floating Body: "<<Mass_fb<<endl;
			cout<<"Density Floating Body: "<<Rfb<<endl;
			cout<<"Moments of Inertia Tensor:\n"<<I_<<endl;
		}
		
}
