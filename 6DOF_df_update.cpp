/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_df.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void sixdof_df::updateFSI(lexer *p, fdm *a, ghostcell* pgc, bool& converged)
{
    // Print quaternion
	if(p->mpirank==0 && converged == true)
    {
		//cout<<"Quaternion: "<<e_(0)<<" "<<e_(1)<<" "<<e_(2)<<" "<<e_(3)<<endl;
    }
	
    // Update transformation matrix (Shivarama PhD thesis, p. 19)
    quat_matrices(e_);


    // Calculate new position
    updatePosition(p, a, pgc, converged);


    // Update angular velocities 
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
      
    
    if (converged == true)
    {
        // Global body variables
        interface(p,false);    
        maxvel(p,a,pgc);
    }
}


void sixdof_df::updatePosition(lexer *p, fdm *a, ghostcell *pgc, bool converged)
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

	if(p->mpirank==0 && converged == true)
    {
        cout<<"XG: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
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
        
            tri_x[n][q] = point(0) + c_(0);
            tri_y[n][q] = point(1) + c_(1);
            tri_z[n][q] = point(2) + c_(2);

			// 2D
			if(p->X11_v != 1 && p->X11_p != 1 && p->X11_r != 1) 
			{
				tri_y[n][q] = tri_y0[n][q] + c_(1);	
			}
        }

	}
	
    // Update floating field
	ray_cast(p,a,pgc);
	reini_AB2(p,a,pgc,a->fb);
    pgc->start4a(p,a->fb,50);    
    
    ULOOP
    {
        a->fbh1(i,j,k) = Hsolidface(p,a,1,0,0);
    }
    VLOOP
    {
        a->fbh2(i,j,k) = Hsolidface(p,a,0,1,0);
    }
    WLOOP
    {
        a->fbh3(i,j,k) = Hsolidface(p,a,0,0,1);
    }
    LOOP
    {
        a->fbh4(i,j,k) = Hsolidface(p,a,0,0,0);
    }

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}


void sixdof_df::finalise(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
    // Store RBM motion
    saveTimeStep(p,alpha);
    
    
    // Store body velocities and position
    interface(p,true);
    
    
    // Print data
    print_stl(p,a,pgc);
    print_parameter(p, a, pgc);
    
    nCorr = 1;
    
    if (p->mpirank == 0)
	{
		cout<<"Ue: "<<p->ufbi<<" Ve: "<<p->vfbi<<" We: "<<p->wfbi<<" Pe: "<<p->pfbi<<" Qe: "<<p->qfbi<<" Re: "<<p->rfbi<<endl;
    }
}


void sixdof_df::quat_matrices(const Eigen::Vector4d& e)
{
    E_ << -e(1), e(0), -e(3), e(2),
         -e(2), e(3), e(0), -e(1),
         -e(3), -e(2), e(1), e(0); 

    G_ << -e(1), e(0), e(3), -e(2),
         -e(2), -e(3), e(0), e(1),
         -e(3), e(2), -e(1), e(0); 

    R_ = E_*G_.transpose(); 
    Rinv_ = R_.inverse();

    p->quatRotMat = R_;
}


double sixdof_df::Hsolidface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_fb;
	
    psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->knoy == 1)
    {
        psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    }

    // Construct solid heaviside function

    phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc));
    
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
