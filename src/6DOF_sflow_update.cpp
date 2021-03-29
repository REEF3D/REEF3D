/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sixdof_sflow::updateFSI(lexer *p, ghostcell* pgc)
{
    // Update transformation matrix (Shivarama PhD thesis, p. 19)
    quat_matrices(e_);

    // Calculate new position
    updatePosition(p, pgc);
      
    // Global body variables
    //interface(p,false);    
    //maxvel(p,a,pgc);
}


void sixdof_sflow::updatePosition(lexer *p, ghostcell *pgc)
{
	// Calculate Euler angles from quaternion
	
	// around z-axis
	psi = atan2(2.0*(e_(1)*e_(2) + e_(3)*e_(0)), 1.0 - 2.0*(e_(2)*e_(2) + e_(3)*e_(3))); 
	
	// around new y-axis
	double arg = 2.0*(e_(0)*e_(2) - e_(1)*e_(3));
	
	if (fabs(arg) >= 1.0)
	{
		theta = SIGN(arg)*PI/2.0;
	}
	else
	{
		theta = asin(arg);														
	}	
		
	// around new x-axis
	phi = atan2(2.0*(e_(2)*e_(3) + e_(1)*e_(0)), 1.0 - 2.0*(e_(1)*e_(1) + e_(2)*e_(2)));

	if(p->mpirank==0)
    {
        cout<<"XG: "<<p->xg<<" YG: "<<p->yg<<" ZG: "<<p->zg<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;

		//cout<<"Ue: "<<p_(0)/Mass_fb<<" Ve: "<<p_(1)/Mass_fb<<" We: "<<p_(2)/Mass_fb<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }

	// Update position of triangles 
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {
            // Update coordinates of triangles 
            // (tri_x0 is vector between tri_x and xg)
  
            Eigen::Vector3d point(tri_x0[n][q], tri_y0[n][q], tri_z0[n][q]);
					
            point = R_*point;
        
            tri_x[n][q] = point(0) + p->xg;
            tri_y[n][q] = point(1) + p->yg;
            tri_z[n][q] = point(2) + p->zg;
        }

	}
	
    // Update floating level set function
	ray_cast(p,pgc);
	reini(p,pgc,fb);
}


void sixdof_sflow::updateForcing_hemisphere(lexer *p, ghostcell *pgc)
{
    // Calculate hemisphere forcing fields
    double H, press0, r, xpos, ypos, dist;

    press0 = 300.0;
    r = 40.0;

	SLICELOOP1
    {
        xpos = p->pos1_x() - p->xg;
        ypos = p->pos1_y() - p->yg;
        dist = xpos*xpos + ypos*ypos;

        H = Hsolidface(p,1,0);
       
        if (dist < r*r)
        {
            press_x(i,j) = -xpos*press0/(r*r*sqrt(1.0 - dist/(r*r))*p->W1);
        }
        else
        {
            press_x(i,j) = 0.0;
        }
    }
    
	SLICELOOP2
    {
        xpos = p->pos2_x() - p->xg;
        ypos = p->pos2_y() - p->yg;
        dist = xpos*xpos + ypos*ypos;

        H = Hsolidface(p,0,1);
       
        if (dist < r*r)
        {
            press_y(i,j) = -ypos*press0/(r*r*sqrt(1.0 - dist/(r*r))*p->W1);
        }
        else
        {
            press_y(i,j) = 0.0;
        }
    }
	
    pgc->gcsl_start1(p,press_x,10);
    pgc->gcsl_start2(p,press_y,11);
};


void sixdof_sflow::updateForcing_ship(lexer *p, ghostcell *pgc)
{
    // Calculate hemisphere forcing fields
    double H, press0, xpos, ypos, Ls, Bs, as, cl, cb;

    press0 = 30.0;
    as = 16.0; 
    cl = 2.0;
    cb = 16.0;

    Ls = 100.0;
    Bs = 20.0;

	SLICELOOP1
    {
        xpos = p->pos1_x();
        ypos = p->pos1_y();
        H = Hsolidface(p,1,0);
        
        press_x(i,j) = -H*4.0*press0*cl/Ls*pow(xpos/Ls,3)*(1.0 - cb*pow(ypos/Bs,2))*exp(-as*ypos*ypos/(Bs*Bs));
    }
    
	SLICELOOP2
    {
        xpos = p->pos2_x();
        ypos = p->pos2_y();
        H = Hsolidface(p,0,1);
        
        press_y(i,j) = -H*2.0*press0*as/Bs*(1.0 - cl*pow(xpos/Ls,4))*(cb/as + 1.0 - cb*pow(ypos/Bs,2))*ypos/Bs*exp(-as*ypos*ypos/(Bs*Bs));
    }
	
    pgc->gcsl_start1(p,press_x,10);
    pgc->gcsl_start2(p,press_y,11);
};


double sixdof_sflow::Hsolidface(lexer *p, int aa, int bb)
{
    double psi, H, phival_fb;

    psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP]); 

    // Construct solid heaviside function

    phival_fb = 0.5*(fb(i,j) + fb(i+aa,j+bb));
    
    if (-phival_fb > psi)
    {
        H = 1.0;
    }
    else if (-phival_fb < -psi)
    {
        H = 0.0;
    }
    else
    {
        H = 0.5*(1.0 + -phival_fb/psi + (1.0/PI)*sin((PI*-phival_fb)/psi));
    }
        
    return H;
}
