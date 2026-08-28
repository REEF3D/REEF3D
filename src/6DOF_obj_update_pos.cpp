/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>

void sixdof_obj::update_position_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalize)
{
    // Calculate new position
    update_Euler_angles(p,pgc);

    // Update STL mesh
    update_trimesh_3D(p,a,pgc,finalize);
    
    // Update angular velocities 
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
    
    if(p->mpirank==0 && finalize==true)
    {
        cout<<"XG["<<n6DOF<<"]: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
        cout<<"Ue: "<<u_fb(0)<<" Ve: "<< u_fb(1)<<" We: "<< u_fb(2)<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }

}

void sixdof_obj::update_Euler_angles(lexer *p, ghostcell *pgc)
{
	// Calculate Euler angles from quaternion
	
	// around z-axis
	psi = atan2(2.0*(e_(1)*e_(2) + e_(3)*e_(0)), 1.0 - 2.0*(e_(2)*e_(2) + e_(3)*e_(3))); 
	
	// around new y-axis
	double arg = 2.0*(e_(0)*e_(2) - e_(1)*e_(3));
	
	if (fabs(arg) >= 1.0)
	theta = SIGN(arg)*PI/2.0;
    
	else
	theta = asin(arg);														
			
	// around new x-axis
	phi = atan2(2.0*(e_(2)*e_(3) + e_(1)*e_(0)), 1.0 - 2.0*(e_(1)*e_(1) + e_(2)*e_(2)));	
}

void sixdof_obj::update_trimesh_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalize)
{
	// Update position of triangles 
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {
            // Update coordinates of triangles 
            Eigen::Vector3d point(tri_x0[n][q], tri_y0[n][q], tri_z0[n][q]);
					
            point = R_*point;
        
            tri_x[n][q] = point(0) + c_(0);
            tri_y[n][q] = point(1) + c_(1);
            tri_z[n][q] = point(2) + c_(2);
        }
	}

    // Update floating level set function
	ray_cast(p,a,pgc);
	reini_RK2(p,a,pgc,a->fb);
    
    pgc->start4a(p,a->fb,50);   
}

void sixdof_obj::union_fb(lexer *p, fdm *a, field &fb_all)
{
    ALOOP
    fb_all(i,j,k) = MIN(fb_all(i,j,k), a->fb(i,j,k));
}

void sixdof_obj::reini_fb(lexer *p, fdm *a, ghostcell *pgc)
{
    reini_RK2(p,a,pgc,a->fb);
    pgc->start4a(p,a->fb,50);
}

void sixdof_obj::add_solid_heaviside(lexer *p, fdm *a, ghostcell *pgc)
{
    ULOOP
    a->fbh1(i,j,k) = MIN(a->fbh1(i,j,k) + Hsolidface(p,a,1,0,0), 1.0);

    VLOOP
    a->fbh2(i,j,k) = MIN(a->fbh2(i,j,k) + Hsolidface(p,a,0,1,0), 1.0);

    WLOOP
    a->fbh3(i,j,k) = MIN(a->fbh3(i,j,k) + Hsolidface(p,a,0,0,1), 1.0);

    LOOP
    a->fbh4(i,j,k) = MIN(a->fbh4(i,j,k) + Hsolidface(p,a,0,0,0), 1.0);
}

void sixdof_obj::apply_kinematics(lexer *p, const Eigen::Vector3d& c, const Eigen::Vector4d& e,
                                  const Eigen::Vector3d& v, const Eigen::Vector3d& omegaI)
{
    c_ = c;
    e_ = e;
    if(e_.norm()>1.0e-14)
    e_.normalize();

    quat_matrices(p);

    p_ = Mass_fb * v;
    omega_I = omegaI;
    omega_B = Rinv_ * omega_I;
    h_ = I_ * omega_B;
}

void sixdof_obj::get_pose_vel(Eigen::Vector3d& c, Eigen::Vector4d& e, Eigen::Vector3d& v, Eigen::Vector3d& omegaI) const
{
    c = c_;
    e = e_;
    if(Mass_fb>1.0e-14)
    v = p_/Mass_fb;
    else
    v.setZero();
    omegaI = omega_I;
}

void sixdof_obj::get_wrench(Eigen::Vector3d& F, Eigen::Vector3d& M) const
{
    F = Ffb_;
    M = Mfb_;
}

void sixdof_obj::collision_extents(double& dx, double& dy, double& dz) const
{
    double xmin=1.0e20, xmax=-1.0e20;
    double ymin=1.0e20, ymax=-1.0e20;
    double zmin=1.0e20, zmax=-1.0e20;

    for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        xmin = std::min(xmin, tri_x0[n][q]);
        xmax = std::max(xmax, tri_x0[n][q]);
        ymin = std::min(ymin, tri_y0[n][q]);
        ymax = std::max(ymax, tri_y0[n][q]);
        zmin = std::min(zmin, tri_z0[n][q]);
        zmax = std::max(zmax, tri_z0[n][q]);
    }

    dx = std::max(xmax-xmin, 1.0e-4);
    dy = std::max(ymax-ymin, 1.0e-4);
    dz = std::max(zmax-zmin, 1.0e-4);
}

int sixdof_obj::copy_collision_mesh(std::vector<double>& xyz) const
{
    xyz.resize(size_t(std::max(tricount, 0)) * 9);

    for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        const size_t i = size_t(n)*9 + size_t(q)*3;
        xyz[i+0] = tri_x0[n][q];
        xyz[i+1] = tri_y0[n][q];
        xyz[i+2] = tri_z0[n][q];
    }

    return tricount;
}

