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


// REEF3D::CORAL: cross-flow vortex-induced vibration of the dynamic mooring lines (Thorsen's synchronization
// model, as MoorDyn v2): the lift at the inner nodes and the shedding phases, which advance with the nodes in the
// Runge-Kutta sub-steps. In its own file for the same reason as coral_mooring_friction.cpp.

#include"coral/mooring/coral_mooring.h"
#include<cmath>

namespace coral
{

// the lift of inner node i of line L: 1/2 rho d |vp| Cl cos(phi) (q x vp) sl with vp the relative flow normal to
// the line (vrn) and the phase phi of the running stage (vc_->ps); the phase rate 2 pi |vp| / d (cF + dF
// sin(phi_yd - phi)) into vc_->pd. phi_yd: the phase of the cross-flow motion, atan2(-ydd / rms(ydd), yd / rms(yd))
// with yd, ydd the node's velocity and last acceleration along q x vp / |vp|. The mean squares relax towards yd^2,
// ydd^2 with the time constant viv_tau (< 0: five shedding periods d / (cF |vp|)) over the sub-step vh_; upd:
// kept (MoorDyn: 500 of its time steps, the mean squares kept once per step)
Vec3 mooring_system::viv_lift(const line_t &L, int i, const Vec3 &vx, const Vec3 &q, const Vec3 &vrn, double sl,
                              bool upd) const
{
    const mooring_props &p = L.p;
    viv_t &V = *vc_;
    const double vp = vrn.norm();
    if(!(vp > 0.0))
    {
        V.pd[i] = 0.0;
        return Vec3::Zero();
    }
    const Vec3 qv = q.cross(vrn);
    const Vec3 e = qv/vp;
    const double yd = vx.dot(e), ydd = V.acc[i].dot(e);
    const double tau = p.viv_tau > 0.0 ? p.viv_tau : 5.0*p.d/(p.cF*vp);
    const double w = std::min(1.0, vh_/tau);
    const double y2 = V.y2[i] + w*(yd*yd - V.y2[i]);
    const double a2 = V.a2[i] + w*(ydd*ydd - V.a2[i]);
    if(upd)
    {
        V.y2[i] = y2;
        V.a2[i] = a2;
    }
    const double phi_yd = (y2 > 0.0 && a2 > 0.0) ? std::atan2(-ydd/std::sqrt(a2), yd/std::sqrt(y2)) : std::atan2(-ydd, yd);
    const double phi = V.ps[i];
    V.pd[i] = 2.0*M_PI*(p.cF + p.dF*std::sin(phi_yd - phi))*vp/p.d;
    return (0.5*p.rho*p.d*vp*p.Cl*std::cos(phi)*sl)*qv;
}

// the phases of the inner nodes of the lines with VIV for stage k (-1: the sub-step's start)
void mooring_system::viv_stage(int k, double c) const
{
    viv_t &V = *vc_;
    for(const line_t &L : lines_)
    {
        if(!(L.p.Cl > 0.0))
            continue;
        for(int i=L.off+1; i<L.off+L.N; ++i)
            V.ps[i] = k < 0 ? V.phi[i] : V.phi[i] + c*V.d[k][i];
    }
}

void mooring_system::viv_keep(int k) const
{
    viv_t &V = *vc_;
    for(const line_t &L : lines_)
    {
        if(!(L.p.Cl > 0.0))
            continue;
        for(int i=L.off+1; i<L.off+L.N; ++i)
            V.d[k][i] = V.pd[i];
    }
}

// RK4 step of the phases (rates of stages 1-3 kept, stage 4 in pd), mapped to 0 .. 2 pi
void mooring_system::viv_advance(double h) const
{
    viv_t &V = *vc_;
    const double tp = 2.0*M_PI;
    for(const line_t &L : lines_)
    {
        if(!(L.p.Cl > 0.0))
            continue;
        for(int i=L.off+1; i<L.off+L.N; ++i)
        {
            double f = V.phi[i] + h/6.0*(V.d[0][i] + 2.0*V.d[1][i] + 2.0*V.d[2][i] + V.pd[i]);
            f -= tp*std::floor(f/tp);
            V.phi[i] = f;
        }
    }
}

// the start: phases 2 pi j / N over the inner nodes j of each line (MoorDyn), mean squares and accelerations zero;
// a node array that grew (a line failure: a new point) keeps the phases of the lines
void mooring_system::viv_start()
{
    const size_t n = r_.size();
    viv_t &V = vv_;
    const bool fresh = V.phi.empty();
    V.phi.resize(n, 0.0);
    V.y2.resize(n, 0.0);
    V.a2.resize(n, 0.0);
    V.acc.resize(n, Vec3::Zero());
    V.ps.assign(n, 0.0);
    V.pd.assign(n, 0.0);
    for(int k=0; k<3; ++k)
        V.d[k].assign(n, 0.0);
    if(!fresh)
        return;
    for(const line_t &L : lines_)
    {
        if(!(L.p.Cl > 0.0))
            continue;
        for(int j=1; j<L.N; ++j)
            V.phi[L.off+j] = 2.0*M_PI*j/double(L.N);
    }
}

}
