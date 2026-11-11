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

#ifndef SEDIMENT_XFILTER_H_
#define SEDIMENT_XFILTER_H_

#include"slice4.h"
#include"increment.h"

class lexer;
class ghostcell;
class sediment_fdm;

// The bed filter (S 100, S 101) of sediment_exner::filter in steps, for grids that exchange their ghost cells
// with other grids between the steps (mesh refinement and multi-block grids, nhflow_amr / sflow_amr), and for
// the curvilinear grid (G 50).
//
// The smoother is the symmetric product S = (Sy Sx + Sx Sy)/2 of the one-dimensional passes of filter():
// conservative (each pass is), and independent of which index family of a curvilinear grid is i (a block
// turned by 90 deg gives the same result). On a uniform grid without masks Sy Sx = Sx Sy (filter()).
//
//   mask(p, amr):     own (the cells filtered here: sediment cells in the S 77 window, amr > 0.5) and m = own;
//                     m is then exchanged: the ghost cells say whether the neighbour takes part
//   weights(p):       the face weights of the own cells
//   copy(p, f, 0/1):  h = f (0), dh = h - f (1); then h or dh is exchanged
//   sweep1(p, 1/2):   a = Sx(h or dh), b = Sy(h or dh); then a and b are exchanged
//   sweep2(p, f, s):  f = (Sy a + Sx b)/2 (s = 1) or f += (Sy a + Sx b)/2 (s = 0) on the own cells
// The filter:  for S100: copy 0, sweep1 1, sweep2 1;  for S101: copy 1, sweep1 2, sweep2 0
class sediment_xfilter : public increment
{
public:
    sediment_xfilter(lexer*);

    void mask(lexer*, sediment_fdm*, slice*);
    void weights(lexer*);
    void copy(lexer*, slice&, int);
    void sweep1(lexer*, int);
    void sweep2(lexer*, slice&, int);

    // the whole filter on one grid with the ghost cells of pgc (gcsl_start4)
    void run(lexer*, ghostcell*, sediment_fdm*, slice&, int, int);

    slice4 own, m, h, dh, a, b;

private:
    slice4 wxm, wxp, wym, wyp;
    double beta;
    int ydir;
};

#endif
