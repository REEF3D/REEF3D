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
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::transform(lexer *p, fdm *a, ghostcell* pgc, bool finalise)
{
    // Update transformation matrix (Shivarama PhD thesis, p. 19)
    quat_matrices(e_);

    // Calculate new position
    update_Position(p, a, pgc, finalise);

    // Update angular velocities 
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
      
    // Global body variables
    interface(p,false);    
    maxvel(p,a,pgc);
}

void sixdof_df_object::get_trans(lexer *p, fdm *a, ghostcell *pgc, Eigen::Vector3d& dp, Eigen::Vector3d& dc, const Eigen::Vector3d& pp, const Eigen::Vector3d& c)
{
    dp = Ffb_; 
    dc = pp/Mass_fb;

	// Prescribed motions
	prescribedMotion(p,a,pgc,dp,dc);
} 

void sixdof_df_object::get_rot(Eigen::Vector3d& dh, Eigen::Vector4d& de, const Eigen::Vector3d& h, const Eigen::Vector4d& e)
{
    // Update Euler parameter matrices
    quat_matrices(e);

    // RHS of e
    de = 0.5*G_.transpose()*I_.inverse()*h;
    
    // RHS of h
    // Transforming torsion into body fixed system (Shivarama and Schwab)
    Gdot_ << -de(1), de(0), de(3),-de(2),
             -de(2),-de(3), de(0), de(1),
             -de(3), de(2),-de(1), de(0); 
   
    dh = 2.0*Gdot_*G_.transpose()*h + Rinv_*Mfb_;
} 

void sixdof_df_object::quat_matrices(const Eigen::Vector4d& e)
{
    E_ << -e(1), e(0), -e(3), e(2),
          -e(2), e(3), e(0), -e(1),
          -e(3), -e(2), e(1), e(0); 

    G_ << -e(1), e(0), e(3), -e(2),
          -e(2), -e(3), e(0), e(1),
          -e(3), e(2), -e(1), e(0); 

    R_ = E_*G_.transpose(); 
    Rinv_ = R_.inverse();

    quatRotMat = R_;
}


