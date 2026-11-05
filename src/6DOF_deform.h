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

#ifndef SIXDOF_DEFORM_H_
#define SIXDOF_DEFORM_H_

#include<vector>

class sixdof_geometry;

//  Prescribed body deformation of a 6DOF body: undulatory swimming (fish, eel), X 195.
//
//  The deformation is a lateral displacement h of the body-frame hull (tri0, relative to the CoG):
//
//      y_B -> y_B + h(x_B,t)          (3D: lateral axis y_B; 2D, x-z plane: z_B)
//
//      h(x,t) = r(t) L A(s) sin(2 pi (s/lambda - f t)) - c0(t) - c1(t) x
//
//  s = (x_head - x)/L in [0,1] along the body axis (head at the largest body-frame x, the fish
//  swims towards +x_B), L the length of the body, A(s) the amplitude envelope (fraction of L):
//
//      type 1  linear (Carling et al. 1998, Kern & Koumoutsakos 2006):  a (s + 0.03125)/1.03125
//      type 2  carangiform (Videler & Hess 1984, Borazjani & Sotiropoulos 2008):
//              (a/0.1)(0.02 - 0.0825 s + 0.1625 s^2)
//      type 3  anguilliform (Tytell 2004, Borazjani & Sotiropoulos 2009):  a exp(s - 1)
//
//  a: amplitude at the tail (fraction of L), lambda: wave length (fraction of L), f: tail-beat
//  frequency [Hz], r(t): ramp from 0 to 1 over t_ramp (0.5 (1 - cos(pi t/t_ramp))).
//
//  Recoil (X 196 1): c0 + c1 x removes the net lateral momentum and the angular momentum of the
//  deformation about the CoG (least squares with the mass distribution of the body along x), so
//  that the deformation neither moves the CoG nor turns the body: the rigid-body equations of the
//  6DOF body stay valid with the fluid loads alone, and translation and rotation of the body follow
//  from the fluid (self-propulsion, recoil sway and yaw).
//
//  The map (x, y) -> (x, y + h(x)) is a shear: the cross-sections keep their area, and the
//  deformation velocity (0, dh/dt(x)) is divergence free, so the direct forcing of the body
//  needs no source term in the pressure Poisson equation. The midline stretches with the slope of h
//  (by about <h_x^2>/2: 1-5 % for the usual amplitudes).
//
//  Types 11-13: the same envelopes on an inextensible midline (Carling et al. 1998, Kern &
//  Koumoutsakos 2006, Gazzola et al. 2011): y(s,t) = r L A(s) sin(2 pi (s/lambda - f t)) along the
//  arc length s, sin theta = dy/ds, the axial position from the inextensibility (the midline keeps
//  its length L, the body contracts axially by int (1 - cos theta) ds), the cross-sections normal
//  to the midline: the material point (x, l) of the undeformed body (l lateral) goes to
//      M(s) + l (sin theta, cos theta)     (axial, lateral),  M(s) = (x + int_0^s (1 - cos theta), y)
//  Recoil: translation (axial and lateral) and rotation in the plane of the undulation from the
//  momentum and angular momentum of the deformation (volume and the first and second moments of l
//  per bin along x). The velocity at a point of the deformed body (direct forcing) from its
//  material coordinates (Newton iterations on the map). The deformation velocity of these types
//  has the divergence -l (d2 theta/ds dt)/(1 - l dtheta/ds) (the sections turn), which the
//  projection of the direct forcing does not keep.

class sixdof_deform
{
public:

    sixdof_deform(int type, double a, double lambda, double f, double tramp, int lateral, bool recoil);

    // body length, head position and mass distribution along x of the undeformed body frame hull
    void setup(const sixdof_geometry&);

    // kinematics at time t: ramp and recoil terms
    void set_time(double t);
    double time() const {return t_;}

    // lateral displacement and lateral velocity at the body-frame axial coordinate x (relative to the
    // CoG), types 1-3
    double disp(double x) const;
    double vel(double x) const;
    
    // all types: displacement and velocity (axial, lateral) of the material point (x, l) of the
    // undeformed body frame, and the velocity at the point (x, l) of the deformed body
    void displacement(double x, double l, double &dx, double &dl) const;
    void velocity(double x, double l, double &vx, double &vl) const;
    void velocity_at(double x, double l, double &vx, double &vl) const;
    bool inextensible() const {return inext_;}
    
    // discrete recoil correction (X 198): rigid motion in the plane of the undulation after the
    // recoil (translation axial/lateral, rotation psi from the axis towards the lateral axis, about
    // the CoG) and its rates
    struct correction {double tx = 0.0, tl = 0.0, psi = 0.0, vx = 0.0, vl = 0.0, psit = 0.0;};
    const correction &corr() const {return cr_;}
    void set_corr(const correction &c) {cr_ = c; corr_on_ = true;}

    // lateral axis of the body frame: 1 y, 2 z
    int lateral() const {return lat_;}

    // diagnostics
    double length() const {return L_;}
    double xhead() const {return xhead_;}
    double xtail() const {return xtail_;}
    double frequency() const {return f_;}
    double tail_amplitude() const {return L_*envelope(1.0);}     // peak-to-peak / 2 at the tail tip (r = 1)
    double ramp() const {return r_;}
    double recoil0() const {return c0_;}
    double recoil1() const {return c1_;}
    double volume() const {return vol_;}

private:

    void displacement_raw(double, double, double&, double&) const;
    void velocity_raw(double, double, double&, double&) const;
    void velocity_at_raw(double, double, double&, double&) const;
    void correct_velocity(double, double, double&, double&) const;
    correction cr_;
    bool corr_on_ = false;
    
    double envelope(double) const;
    double envelope_s(double) const;
    void raw(double, double&, double&) const;       // h and dh/dt without the recoil terms
    void fit(const std::vector<double>&, double&, double&) const;

    const int type_;
    const bool inext_;
    const int lat_;
    const double a_, lambda_, f_, tramp_;
    const bool recoil_;

    double t_ = 0.0, r_ = 0.0, rdot_ = 0.0;
    double c0_ = 0.0, c1_ = 0.0, d0_ = 0.0, d1_ = 0.0;

    double xhead_ = 0.0, xtail_ = 0.0, L_ = 1.0, vol_ = 0.0;

    // mass (volume) distribution along x: bin centres and volumes
    std::vector<double> xb_, mb_;
    double S0_ = 0.0, S1_ = 0.0, S2_ = 0.0;
    
    // ... and the first and second moments of the lateral coordinate per bin, the moment of inertia
    // of the volume in the plane of the undulation
    std::vector<double> n1b_, n2b_;
    double In_ = 0.0;
    
    // inextensible midline (types 11-13): tables over the arc length, recoil (translation cx, cl,
    // rotation phi and their rates; phi integrated in time)
    struct station {double th, tht, ths, y, yt, ax, axt;};
    struct mtab
    {
        std::vector<double> th, tht, ym, ymt, ax, axt;
        double cx, cl, vx, vl, phit;
    };
    void midline(double, double, double, mtab&) const;
    void ramp(double, double&, double&) const;
    void at(double, station&) const;
    void at(const mtab&, double, station&) const;
    double phit_at(double) const;
    double phi_at(double);
    mtab tab_;
    std::vector<double> phik_;      // recoil rotation at the multiples of T/64
    double cx_ = 0.0, cl_ = 0.0, phi_ = 0.0, vx_ = 0.0, vl_ = 0.0, phit_ = 0.0;
};

#endif
