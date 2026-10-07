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

#ifndef BACKGROUND_STATE_H_
#define BACKGROUND_STATE_H_

#include<vector>

class lexer;

/*--------------------------------------------------------------------
background_state: tidal / current background of the boundary zones
(iowave redesign, step 4).

  B 510 id mode dir t_ramp     mode 1: harmonic, 2: time series, 3: constant;
                               dir [deg]: direction of propagation; t_ramp [s]: spin-up
  B 511 id a T phase           harmonic constituent (repeatable): a cos(2 pi t/T - phase)
  B 514 id eta0 U V            constant level offset [m] and current [m/s]
  B 515 id x0 y0 h_ref         progressive background: the tide travels in direction dir
                               as a long wave on the depth h_ref (<= 0: still water level
                               F 60) from (x0,y0); without B 515 the background is
                               uniform in space
  file background-<id>.dat     mode 2: lines "t eta" or "t eta U V", linear in time
  B 530 mode N                 waves on the background (iowave, linear waves): every N steps
                               k of each source from the background of its zones,
                               mode 1: on h_eff = h + eta_b, 2: h_eff and Doppler
                               omega = sigma + k U_n (iowave_nhflow_tide.cpp)

  eta_b = r(t) eta0 + eta_t,   eta_t = r(t) sum a cos(2 pi t/T - phase - k s)  (mode 1)
                                     = r(t) eta_file(t - s/c)                  (mode 2)
  U_b   = r(t) U + eta_t sqrt(g/h) cos(dir)   (or r(t) U_file(t - s/c) if the file has U, V)
  V_b   = r(t) V + eta_t sqrt(g/h) sin(dir)

with s = (x-x0) cos(dir) + (y-y0) sin(dir) (0 without B 515), c = sqrt(g h_ref),
k = 2 pi/(T c), h the local still water depth and r(t) the spin-up ramp.
--------------------------------------------------------------------*/

class background_state
{
public:
    void read(lexer*);
    void update(lexer*, double);            // evaluates the time functions at time t
    
    bool empty() const {return bg.empty();}
    int index(int) const;                   // of a background id; -1: none
    
    double eta(int b, double x, double y) const;
    void vel(int b, double h, double x, double y, double &u, double &v) const;   // depth averaged
    
private:
    struct constituent {double a, T, phase, k, C=1.0, S=0.0;};
    struct item
    {
        int id=0, mode=0;
        double dir=0.0, tramp=0.0;
        double eta0=0.0, U=0.0, V=0.0;
        std::vector<constituent> c;
        bool prog=false, file_uv=false;
        double x0=0.0, y0=0.0, href=0.0, cel=0.0;
        std::vector<double> ft, feta, fu, fv;   // time series (mode 2)
        double r=1.0;                            // at the time of the last update
    };
    std::vector<item> bg;
    double g=9.81, time=0.0;
    
    double tide(const item&, double, double, double&, double&) const;   // eta_t, file U, V
    static double interp(const std::vector<double>&, const std::vector<double>&, double);
};

#endif
