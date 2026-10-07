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

void sixdof_actuator_disk::accumulate(const Eigen::Vector3d &x, double V, double fl, int c, sums &s) const
{
    double wa, wt, r;
    Eigen::Vector3d et;
    
    if(!weights(x,wa,wt,et,r))
    return;
    
    // lever of a force along the unit vector of component c about the axis
    const double g = (x - centre).cross(Eigen::Vector3d::Unit(c)).dot(axis);
    
    s.SA0 += wa*V;
    s.SA  += fl*wa*V;
    s.SW  += fl*wt*V;
    s.SE  += fl*wt*et(c)*V;
    s.G1  += fl*wt*g*et(c)*V;
    s.G0  += fl*wt*g*V;
    s.GA  += fl*wa*g*V;
}

double sixdof_actuator_disk::swirl_factor(const sums *s) const
{
    // tangential force density kappa fl wt (et_c - SE_c/SW_c): no net force in component c, and
    // the torque about the axis of all components together, thrust part included, is sense Q
    double den = 0.0, tT = 0.0;
    
    for(int c=0; c<3; ++c)
    {
        if(s[c].SW>0.0)
        den += s[c].G1 - s[c].SE/s[c].SW*s[c].G0;
        
        if(s[c].SA>0.0)
        tT -= T*axis(c)/s[c].SA*s[c].GA;
    }
    
    return fabs(den)>1.0e-30 ? (double(sense)*Q - tT)/den : 0.0;
}

double sixdof_actuator_disk::force(const Eigen::Vector3d &x, double fl, int c, const sums &s, double kappa) const
{
    double wa, wt, r;
    Eigen::Vector3d et;
    
    if(!weights(x,wa,wt,et,r))
    return 0.0;
    
    double f = 0.0;
    
    if(s.SA>0.0)
    f -= fl*T*wa/s.SA*axis(c);
    
    if(s.SW>0.0)
    f += kappa*fl*wt*(et(c) - s.SE/s.SW);
    
    return f;
}
