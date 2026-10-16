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

#include"coral/core/coral_rotation.h"
#include<cmath>

namespace coral
{

Mat3 skew(const Vec3 &a)
{
    Mat3 S;
    S <<  0.0, -a(2),  a(1),
         a(2),   0.0, -a(0),
        -a(1),  a(0),   0.0;
    return S;
}

Mat3 exp_so3(const Vec3 &psi)
{
    const double phi2 = psi.squaredNorm();
    const Mat3 S = skew(psi);
    double a, b;    // R = I + a S + b S^2
    if(phi2 < 1.0e-8)
    {
        a = 1.0 - phi2/6.0 + phi2*phi2/120.0;
        b = 0.5 - phi2/24.0 + phi2*phi2/720.0;
    }
    else
    {
        const double phi = std::sqrt(phi2);
        a = std::sin(phi)/phi;
        b = (1.0 - std::cos(phi))/phi2;
    }
    return Mat3::Identity() + a*S + b*S*S;
}

Vec3 log_so3(const Mat3 &R)
{
    // through the quaternion: robust over the whole range [0, pi]
    Eigen::Quaterniond q(R);
    q.normalize();
    if(q.w() < 0.0)
        q.coeffs() *= -1.0;
    const Vec3 v = q.vec();
    const double s = v.norm();
    if(s < 1.0e-12)
        return 2.0*v;    // small angle: psi = 2 v (1 + O(s^2))
    const double phi = 2.0*std::atan2(s, q.w());
    return (phi/s)*v;
}

Mat3 tangent_so3(const Vec3 &psi)
{
    const double phi2 = psi.squaredNorm();
    const Mat3 S = skew(psi);
    double a, b;    // T = I - a S + b S^2
    if(phi2 < 1.0e-8)
    {
        a = 0.5 - phi2/24.0 + phi2*phi2/720.0;
        b = 1.0/6.0 - phi2/120.0 + phi2*phi2/5040.0;
    }
    else
    {
        const double phi = std::sqrt(phi2);
        a = (1.0 - std::cos(phi))/phi2;
        b = (phi - std::sin(phi))/(phi2*phi);
    }
    return Mat3::Identity() - a*S + b*S*S;
}

Vec3 euler_from_matrix(const Mat3 &R)
{
    Vec3 rpy;
    rpy(1) = std::asin(std::max(-1.0, std::min(1.0, -R(2,0))));
    rpy(0) = std::atan2(R(2,1), R(2,2));
    rpy(2) = std::atan2(R(1,0), R(0,0));
    return rpy;
}

Mat3 matrix_from_euler(const Vec3 &rpy)
{
    return (Eigen::AngleAxisd(rpy(2), Vec3::UnitZ())
          * Eigen::AngleAxisd(rpy(1), Vec3::UnitY())
          * Eigen::AngleAxisd(rpy(0), Vec3::UnitX())).toRotationMatrix();
}

Mat3 orthonormalize(const Mat3 &R)
{
    Eigen::Quaterniond q(R);
    q.normalize();
    return q.toRotationMatrix();
}

}
