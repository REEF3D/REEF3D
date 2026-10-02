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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_pto_elements.h"

// --- linear damper ---

pto_damper::pto_damper(double B_) : B(B_)
{
}

pto_output pto_damper::force(const pto_state &s)
{
    pto_output o;

    o.F     = -B*s.qd;
    o.dFdqd = -B;
    o.P     = -o.F*s.qd;

    return o;
}

void pto_damper::set_param(const std::string &key, double val)
{
    if(key=="B")
    B = val;
}

// --- spring-damper ---

pto_springdamper::pto_springdamper(double K_, double B_) : K(K_), B(B_)
{
}

pto_output pto_springdamper::force(const pto_state &s)
{
    pto_output o;

    o.F     = -K*s.q - B*s.qd;
    o.dFdq  = -K;
    o.dFdqd = -B;
    o.P     = -o.F*s.qd;

    return o;
}

void pto_springdamper::set_param(const std::string &key, double val)
{
    if(key=="K")
    K = val;

    if(key=="B")
    B = val;
}

// --- end-stops ---

pto_endstop::pto_endstop(double qmin_, double qmax_, double K_, double C_) : qmin(qmin_), qmax(qmax_), K(K_), C(C_)
{
}

pto_output pto_endstop::force(const pto_state &s)
{
    // penetration d beyond the stroke; the damping part may not pull the body back
    // into the stop when it moves out of it (F is clipped to zero then)
    pto_output o;

    if(s.q>qmax)
    {
        const double F = -K*(s.q - qmax) - C*s.qd;

        if(F<0.0)
        {
        o.F     = F;
        o.dFdq  = -K;
        o.dFdqd = -C;
        }
    }

    else if(s.q<qmin)
    {
        const double F = -K*(s.q - qmin) - C*s.qd;

        if(F>0.0)
        {
        o.F     = F;
        o.dFdq  = -K;
        o.dFdqd = -C;
        }
    }

    o.P = -o.F*s.qd;

    return o;
}
