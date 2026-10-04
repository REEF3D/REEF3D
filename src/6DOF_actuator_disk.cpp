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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_load.h"
#include"definitions.h"
#include<cmath>

bool sixdof_actuator_disk::weights(const Eigen::Vector3d &x, double &wa, double &wt, Eigen::Vector3d &et, double &r) const
{
    wa = wt = r = 0.0;
    
    const Eigen::Vector3d d = x - centre;
    const double a = d.dot(axis);
    
    if(fabs(a) > 0.5*thickness || R<=0.0)
    return false;
    
    const Eigen::Vector3d rv = d - a*axis;
    r = rv.norm();
    
    if(r>=R || r<=Rh)
    return false;
    
    // Hough-Ordway radial distribution (Stern et al. 1988), r* in [0,1] from hub to tip
    const double rh = Rh/R;
    const double rs = (r/R - rh)/(1.0 - rh);
    const double shape = rs*sqrt(1.0 - rs);
    
    // smooth axial profile over the thickness of the source region
    const double ca = cos(PI*a/thickness);
    const double fa = ca*ca;
    
    wa = shape*fa;
    wt = shape/(rs*(1.0 - rh) + rh)*fa;
    
    et = double(sense)*axis.cross(rv/r);
    
    return true;
}
