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
#include"fdm.h"
#include"ghostcell.h"
#include"6DOF_motionext.h"

void sixdof_obj::get_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d& dp_, Eigen::Vector3d& dc_, Eigen::Vector3d& pp_, Eigen::Vector3d& c_)
{
    dp_ = Ffb_; 
    dc_ = pp_/Mass_fb;

	// External motions
	pmotion->motionext_trans(p,pgc,dp_,dc_);
} 

void sixdof_obj::get_rot(lexer *p, Eigen::Vector3d& dh, Eigen::Vector4d& de_, Eigen::Vector3d& h_, Eigen::Vector4d& e_)
{
    // Update Euler parameter matrices
    quat_matrices();

    // RHS of e
    de_ = 0.5*G_.transpose()*I_.inverse()*h_;
    
    // RHS of h
    // Transforming torsion into body fixed system (Shivarama and Schwab)
    Gdot_ << -de_(1), de_(0), de_(3),-de_(2),
             -de_(2),-de_(3), de_(0), de_(1),
             -de_(3), de_(2),-de_(1), de_(0); 
   
    dh_ = 2.0*Gdot_*G_.transpose()*h_ + Rinv_*Mfb_;
    
    // External motions
    pmotion->motionext_rot(p,dh_,h_,de_,G_,I_);
} 

void sixdof_obj::quat_matrices()
{
    // Update transformation matrix (Shivarama PhD thesis, p. 19)
    E_ << -e_(1), e_(0), -e_(3), e_(2),
          -e_(2), e_(3), e_(0), -e_(1),
          -e_(3), -e_(2), e_(1), e_(0); 

    G_ << -e_(1), e_(0), e_(3), -e_(2),
          -e_(2), -e_(3), e_(0), e_(1),
          -e_(3), e_(2), -e_(1), e_(0); 

    R_ = E_*G_.transpose(); 
    Rinv_ = R_.inverse();

    quatRotMat = R_;
}


