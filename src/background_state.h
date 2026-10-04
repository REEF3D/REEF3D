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
(iowave redesign, step 4). A background is uniform in space; zones on
different edges use different backgrounds (e.g. the phase of the tide at
each end of a channel).

  B 510 id mode dir t_ramp     mode 1: harmonic, 3: constant; dir [deg]: direction of
                               propagation of the harmonic tide; t_ramp [s]: spin-up
  B 511 id a T phase           harmonic constituent (repeatable): a cos(2 pi t/T - phase)
  B 514 id eta0 U V            constant level offset [m] and current [m/s]

  eta_b = r(t) (eta0 + sum a cos(2 pi t/T - phase))
  U_b   = r(t) U + eta_h sqrt(g/h) cos(dir),  V_b = r(t) V + eta_h sqrt(g/h) sin(dir)

with eta_h the harmonic part (progressive long wave, depth-uniform),
h the still water depth and r(t) the spin-up ramp.
--------------------------------------------------------------------*/

class background_state
{
public:
    void read(lexer*);
    void update(lexer*, double);            // evaluates all backgrounds at time t
    
    bool empty() const {return bg.empty();}
    int index(int) const;                   // of a background id; -1: none
    
    double eta(int b) const {return bg[b].eta;}
    void vel(int b, double h, double &u, double &v) const;   // depth averaged, still water depth h
    
private:
    struct constituent {double a, T, phase;};
    struct item
    {
        int id=0, mode=0;
        double dir=0.0, tramp=0.0;
        double eta0=0.0, U=0.0, V=0.0;
        std::vector<constituent> c;
        double eta=0.0, eta_h=0.0, r=1.0;   // at the time of the last update
    };
    std::vector<item> bg;
    double g=9.81;
};

#endif
