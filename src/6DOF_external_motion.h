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
Author: Hans Bihs
--------------------------------------------------------------------*/

class lexer;
class fdm;
class fdm_nhf;
class fdm2D;
class ghostcell;
class vrans;
class net;
class field;
#include <Eigen/Dense>

using namespace std;

#ifndef SIXDOF_EXTERNAL_MOTION_H_
#define SIXDOF_EXTERNAL_MOTION_H_

class sixdof_external_motion
{
public:
    virtual void external_motion_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&)=0;
    virtual void external_motion_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&)=0;

    virtual void ini(lexer*,ghostcell*)=0;
};

#endif
