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

#include"6DOF_sflow.h"
#include"lexer.h"
#include"ghostcell.h"

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
    // Calculate hemisphere pressure field
    double H, press0, r, xpos, ypos, dist;

    press0 = p->X401_p0;
    r = p->X133_rad;
   
    SLICELOOP4
    {
        xpos = p->pos_x() - p->xg;
        ypos = p->pos_y() - p->yg;
        dist = xpos*xpos + ypos*ypos;
        H = Hsolidface(p,0,0);
       
        if (dist < r*r)
        {
            press(i,j) = -H*press0*sqrt(1.0 - dist/(r*r))*ramp_draft(p);
        }
        else
        {
            press(i,j) = 0.0;
        }
    }
    
    pgc->gcsl_start4(p,press,50);
}

void sixdof_sflow::updateForcing_box(lexer *p, ghostcell *pgc)
{
    // Calculate ship-like pressure field
    double H, press0, xpos, ypos, Ls, Bs, as, cl, cb;

    press0 = p->X401_p0;
    as = p->X401_a; 
    cl = p->X401_cl;
    cb = p->X401_cb;

   
    Ls = p->X110_xe[0] - p->X110_xs[0];
    Bs = p->X110_ye[0] - p->X110_ys[0];

	SLICELOOP4
    {
        xpos = p->pos_x() - p->xg;
        ypos = p->pos_y() - p->yg;
        H = Hsolidface(p,0,0);
        
        if (xpos <= Ls/2.0 && xpos >= -Ls/2.0 && ypos <= Bs/2.0 && ypos >= -Bs/2.0)
        {
            press(i,j) = -H*press0*(1.0 - cl*pow(xpos/Ls,4.0))*(1.0 - cb*pow(ypos/Bs,2.0))*exp(-as*pow(ypos/Bs,2.0))*ramp_draft(p);
        }
        else
        {
            press(i,j) = 0.0;
        }
    }
    
    pgc->gcsl_start4(p,press,50);
}

void sixdof_sflow::updateForcing_stl(lexer *p, ghostcell *pgc)
{
    // Calculate ship-like pressure field
    double H, press0, xpos, ypos, as, cl, cb;

    press0 = p->X401_p0;
    as = p->X401_a; 
    cl = p->X401_cl;
    cb = p->X401_cb;

	SLICELOOP4
    {
        H = Hsolidface(p,0,0);

        press(i,j) = -H*fabs(p->W22)*p->W1*draft(i,j)*ramp_draft(p);
    }
    
    pgc->gcsl_start4(p,press,50);
}

void sixdof_sflow::updateForcing_oned(lexer *p, ghostcell *pgc)
{
    // Calculate 1D pressure field
    double press0, xpos, as;

    press0 = p->X401_p0;
    as = p->X401_a; 

	SLICELOOP4
    {
        xpos = p->pos_x() - p->xg;
   
        press(i,j) = press0*exp(-pow(xpos/as,2));
    }

    pgc->gcsl_start4(p,press,50);
}

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
