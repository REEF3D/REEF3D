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

void sixdof_obj::geometry_f(double& w0, double& w1, double& w2, double& f1, double& f2, double& f3, double& g0, double& g1, double& g2)
{
	double temp0 = w0 + w1;
    f1 = temp0 + w2;
    double temp1 = w0 * w0;
    double temp2 = temp1 + w1 * temp0;
    f2 = temp2 + w2 * f1;
    f3 = w0 * temp1 + w1 * temp2 + w2 * f2; 
    g0 = f2 + w0 * (f1 + w0);
    g1 = f2 + w1 * (f1 + w1);
    g2 = f2 + w2 * (f1 + w2);
}

void sixdof_obj::geometry_stl(lexer *p, ghostcell *pgc)
{
    if(p->X180==1)
    {
        double x1, x2, x3, y1, y2, y3, z1, z2, z3;
        Vfb=0.0;
        for (int n = 0; n < tricount; ++n)
        {
        
            
            x1 = tri_x[n][0];
            x2 = tri_x[n][1];
            x3 = tri_x[n][2];
            
            y1 = tri_y[n][0];
            y2 = tri_y[n][1];
            y3 = tri_y[n][2];
            
            z1 = tri_z[n][0];
            z2 = tri_z[n][1];
            z3 = tri_z[n][2];  
                
            Vfb += (1.0/6.0)*(-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
        }
        
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
}

void sixdof_obj::geometry_ls(lexer *p, fdm *a, ghostcell *pgc)
{
	// Total Volume
	double H;
	Vfb=0.0;
   
	ALOOP
	{
        H = Hsolidface(p,a,0,0,0);
		Vfb+= p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H;
	}
	Vfb=pgc->globalsum(Vfb);
   
    // Mass and density calculation
	if(p->X21==1)
	{
	Rfb = p->X21_d;
	Mass_fb = Vfb*Rfb;
    p->X22_m = Mass_fb;
	}
	
	if(p->X22==1)
	{
	Mass_fb = p->X22_m;
	Rfb = Mass_fb/Vfb;
	}
	
	if(p->mpirank==0)
	{
	cout<<"Volume Floating Body: "<<Vfb<<endl;
	cout<<"Mass Floating Body: "<<Mass_fb<<endl;
	cout<<"Density Floating Body: "<<Rfb<<endl;
	}

    p->W_fb = Rfb;


	
	// Origin
	
	double xorig=0.0;
	double yorig=0.0;
	double zorig=0.0;
	double rx, ry, rz, Vol;

	// Center of Gravity
	
	if (p->X23==0)
	{
		
		c_(0) = c_(1) = c_(2) = 0.0;
		
		LOOP
		{
            rx=p->pos_x()-xorig;
            ry=p->pos_y()-yorig;
            rz=p->pos_z()-zorig;
            
            //cout<<p->mpirank<<" CG: "<<rx<<" "<<ry<<" "<<rz<<" . "<<a->fb(i,j,k)<<endl;
            
            H = Hsolidface(p,a,0,0,0);

		    Vol = H*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

            c_(0) += (1.0/Mass_fb)*rx*Rfb*Vol;
            c_(1) += (1.0/Mass_fb)*ry*Rfb*Vol;
            c_(2) += (1.0/Mass_fb)*rz*Rfb*Vol;
		}
	
        c_(0) = pgc->globalsum(c_(0));
        c_(1) = pgc->globalsum(c_(1));
        c_(2) = pgc->globalsum(c_(2));
	
        c_(0) += xorig;
        c_(1) += yorig;
        c_(2) += zorig;
	}
    else if (p->X23==1)
	{
		c_(0) = p->X23_x; 
		c_(1) = p->X23_y; 
		c_(2) = p->X23_z; 
	}
	
	if(p->mpirank==0)
	cout<<"Center of Gravity   xg: "<<c_(0)<<" yg: "<<c_(1)<<" zg: "<<c_(2)<<endl;


	double xgn = c_(0);
	double ygn = c_(1);
	double zgn = c_(2);


// Moments of Inertia

    double Ix,Iy,Iz;	
        
	if(p->X24==0)
	{
        Ix = Iy = Iz = 0.0;
		
		
		// Non-homogenious body just needs rho as a function of space, integration is the same
		
		LOOP
		{
			rx=p->pos_x() - c_(0);
			ry=p->pos_y() - c_(1);
			rz=p->pos_z() - c_(2);
			
            H = Hsolidface(p,a,0,0,0);
			
            Vol = H*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

			Ix += (ry*ry + rz*rz)*Rfb*Vol;
			Iy += (rx*rx + rz*rz)*Rfb*Vol;
			Iz += (rx*rx + ry*ry)*Rfb*Vol;
		
			I_(0,1) -= (rx*ry)*Rfb*Vol;
			I_(0,2) -= (rx*rz)*Rfb*Vol;
			I_(1,2) -= (ry*rz)*Rfb*Vol;
		}

		Ix=pgc->globalsum(Ix);
		Iy=pgc->globalsum(Iy);
		Iz=pgc->globalsum(Iz);

		I_(0,1) = pgc->globalsum(I_(0,1));
		I_(0,2) = pgc->globalsum(I_(0,2));
		I_(1,2) = pgc->globalsum(I_(1,2));	
	}
    else if (p->X24==1)
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
	
	if(p->mpirank==0)
	cout<<"Moments of Inertia Tensor:\n"<<I_<<endl;
}
