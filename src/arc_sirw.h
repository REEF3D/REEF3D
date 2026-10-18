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

#ifndef ARC_SIRW_H_
#define ARC_SIRW_H_

#include"increment.h"
#include<vector>

class lexer;
class ghostcell;
class slice;

// Real-time separation of incident and reflected waves with two gauges and FIR filters (SIRW, Frigaard &
// Brorsen 1995, Aalborg University) for active wave absorption at a boundary (B 98 8 at the wavemaker x-,
// B 99 8 at the beach x+), FNPF and NHFLOW, flume (y-averaged) version.
//
// Two gauges at x1, x2 = x1 + dx in front of the boundary (distance d). Per frequency the complex responses
// H1, H2 that give the wave travelling towards the boundary at the boundary position xb from the two gauge
// signals follow from
//     A_m = I exp(-i k x_m) + R exp(i k x_m),  m = 1,2:   H1 A1 + H2 A2 = R exp(i k xb)   (x-, reflected wave)
//                                                        H1 A1 + H2 A2 = I exp(-i k xb)  (x+, incident wave)
// with |1/(2 sin k dx)| capped at 5 near the singularities and a cosine taper of the band flow..fhigh. The
// horizontal velocity of that progressive wave, u(z) = -+ eta omega cosh(k(z+h))/sinh(kh) (shallow-water limit
// eta sqrt(g/h)), is folded into the responses on nlev levels between bed and still-water level, so the FIR
// filters give the boundary velocity correction per level directly (interpolated in z, the top level above); the
// same for the surface elevation and the linear non-hydrostatic pressure p_nh/rho = g eta (cosh(k(z+h))/cosh(kh) - 1)
// (NHFLOW: Dirichlet value of the Poisson equation at the edge). The filters (frequency sampling, Delta t_f = 1/(8 fhigh)) keep
// only the causal part: the wave passes the gauges before it reaches the boundary, d of about 0.7 L(flow)
// (Frigaard & Brorsen) makes the non-causal part negligible.
// The gauges are sampled once per time step, interpolated to the filter rate, averaged over y (all ranks).
class arc_sirw : public increment
{
public:
    // side -1: x- (wavemaker, reflected wave), +1: x+ (beach, incident wave); xb boundary position, h depth,
    // d distance of the first gauge from the boundary, dx gauge spacing, flow/fhigh band [Hz]
    void ini(lexer*, int side, double xb, double h, double d, double dx, double flow, double fhigh);
    
    // once per time step with the current free surface: sample, filter
    void update(lexer*, ghostcell*, slice &eta, double t);
    
    // boundary velocity correction [m/s] at time t and elevation z above still water (linear extrapolation over
    // at most one filter step, linear interpolation between the levels)
    double velocity(double t, double z) const;
    
    // non-hydrostatic pressure / density [m2/s2] of that wave at elevation z (NHFLOW: Dirichlet value of the
    // Poisson equation at the boundary)
    double pressure(double t, double z) const;
    
    // surface elevation [m] of the wave towards the boundary at the boundary (same filters without the velocity
    // profile), for the shallow-water feedback on the remaining deviation
    double eta(double t) const;
    
    bool on=false;
    double x1=0.0,x2=0.0,flow=0.0,fhigh=0.0,dtf=0.0;
    int M=0;
    
private:
    double gauge(lexer*, ghostcell*, slice&, double);
    
    int side=0;
    int count_last=-1;
    int nlev=11, nch=23;
    double hd=1.0;
    std::vector<std::vector<double> > h1,h2;      // [channel][tap]: velocity per level, surface elevation, pressure per level
    std::vector<double> b1,b2;
    int head=0, nsample=0;
    double tnext=0.0, tprev=-1.0, e1prev=0.0, e2prev=0.0;
    std::vector<double> ycur, yprev;
    double tcur=-1.0;
};

#endif
