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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#ifndef NHFLOW_SHIPBL_H_
#define NHFLOW_SHIPBL_H_

#include<vector>
#include<array>
#include<string>

// Integral boundary layer and displacement body of a ship hull given by its surface triangles.
// C++ version of Cases/KCS/shipbl.py (20 Nov 2026), which stays the reference: the steps and their constants are the
// same, and every step was compared with it (KCS model 1).
//
//   set_hull     triangles of the hull in the tank frame; hull frame: x from the FP to the stern, y to port from the
//                centre plane, z up from the still water level
//   build_grid   surface grid (s, t): stations s = x/L every ds, nz cells of equal arc length along the wetted port
//                half of the section from the centre plane at the depth z_env(x) to the still water line (where
//                there is no hull: the centre plane itself, i.e. the wake)
//   double_body  Hess & Smith source panels of the double model (port unknowns, three images) -> velocity at the wall
//   solve        Thwaites to the trip, then Head + Ludwieg-Tillmann along the outer streamlines with their
//                convergence, cross transport and girthwise diffusion; wake with cf = 0; surface wake = the layer
//                that leaves through the still water line (model assumption, fsw)
//   body         displacement body: every section moved along its normal by delta*, wake behind the stern
//
// The edge velocity is that of the double body alone ("run-free" body).

class nhflow_shipbl
{
public:
    struct P2 {double y,z;};
    typedef std::array<double,9> tri3;      // x0 y0 z0 x1 y1 z1 x2 y2 z2
    typedef std::vector<P2> curve;

    nhflow_shipbl();

    // input
    double L,T,B;               // length between perpendiculars, draft, beam [m]
    double U,nu;                // speed [m/s], kinematic viscosity [m^2/s]
    double x_fp,y_c,z_w;        // forward perpendicular, centre plane, still water level in the tank frame [m]
    double xtr,cwake,cdiff;     // trip at s = xtr; relaxation of the wake shape factor; girthwise diffusion
    bool fsw;                   // surface wake in the displacement body
    double ds,s_end;            // station spacing / L, last station / L
    int nz;                     // cells along the girth
    std::string dumpdir;        // not empty: the arrays of every step as raw doubles (checks against shipbl.py)

    void set_hull(const std::vector<tri3>&);
    void build_grid();
    void double_body();
    void solve();
    void body(std::vector<tri3>&);
    void run(const std::vector<tri3> &hull, std::vector<tri3> &dispbody);

    std::string log;            // summary lines of the steps
    double CF,CV,Sw;            // friction and viscous resistance coefficients of the layer, wetted surface of the grid

private:
    // sections
    std::vector<tri3> tri;      // hull frame
    double xmin,xmax;
    std::vector<double> env_x,env_z;
    double z_env(double) const;
    void parts(double x, std::vector<curve>&) const;
    void wet(double x, curve&, std::vector<char>&) const;
    void full(double x, double z_top, curve&, std::vector<char>&) const;

    // grid
    int ns;
    std::vector<double> s,girth;
    std::vector<P2> Pn,Pc,Hn,Hc;             // (ns, nz+1), (ns, nz), (ns+1, nz+1), (ns+1, nz)
    std::vector<char> wall;                  // (ns, nz)
    std::vector<double> e1y,e1z,e2y,e2z,nx,ny,nzc,J,h2;     // (ns, nz)
    void metrics();

    // edge velocity / U at the cell centres, (ns, nz)
    std::vector<double> Vx,Vy,Vz;

    // boundary layer, (ns, nz) and (ns)
    std::vector<double> th,H,dstar,cf,Ue,dplane,Ms,Ds;

    void dump(const char *name, const std::vector<double>&) const;
    void dump(const char *name, const std::vector<P2>&) const;
    void say(const char *fmt, ...);
};

#endif
