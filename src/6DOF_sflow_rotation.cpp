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

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"

void sixdof_sflow::rotation_tri
(
    lexer *p,
    double phi_,double theta_,double psi_, 
    double &xvec,double &yvec,double &zvec, 
    const double& x0, const double& y0, const double& z0
)
{
	// Distance to origin
    double dx = xvec - x0;
    double dy = yvec - y0;
    double dz = zvec - z0;

	// Rotation using Goldstein page 603 (but there is wrong result)
    xvec = dx*(cos(psi_)*cos(theta_)) + dy*(cos(theta_)*sin(psi_)) - dz*sin(theta_);
    yvec = dx*(cos(psi_)*sin(phi_)*sin(theta_)-cos(phi_)*sin(psi_)) + dy*(cos(phi_)*cos(psi_)+sin(phi_)*sin(psi_)*sin(theta_)) + dz*(cos(theta_)*sin(phi_));
    zvec = dx*(sin(phi_)*sin(psi_)+cos(phi_)*cos(psi_)*sin(theta_)) + dy*(cos(phi_)*sin(psi_)*sin(theta_)-cos(psi_)*sin(phi_)) + dz*(cos(phi_)*cos(theta_));
    
	// Moving back
    xvec += x0;
    yvec += y0;
    zvec += z0;
}	

void sixdof_sflow::quat_matrices(const Eigen::Vector4d& e)
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
