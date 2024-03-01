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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_motionext_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

sixdof_motionext_void::sixdof_motionext_void(lexer *p, ghostcell *pgc)
{
}
    
sixdof_motionext_void::~sixdof_motionext_void()
{
}

void sixdof_motionext_void::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_motionext_void::motionext_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d& dp, Eigen::Vector3d& dc)
{
    

}

void sixdof_motionext_void::motionext_rot(lexer *p, Eigen::Vector3d& dh, Eigen::Vector3d& h, Eigen::Vector4d& de, Eigen::Matrix<double, 3, 4>&G_,  Eigen::Matrix3d&I_)
{

}

