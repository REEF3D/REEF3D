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

#ifndef BC_ZONE_H_
#define BC_ZONE_H_

#include<vector>

class lexer;
class ghostcell;

/*--------------------------------------------------------------------
bc_zone: one boundary zone (iowave redesign, phase 2a). Solver-agnostic:
geometry, enforcement method, relaxation weight, distance, and the box
that AMR must not refine.

Geometry (as iowave's B 107 / B 108 zones): a reference line s -> e and a
half width d. The zone is the quadrilateral reaching d to both sides of
the line; the weight depends on the distance to the line.

bc_zone_set holds all zones of a case. from_legacy() translates B 96,
B 98, B 99, B 107 and B 108 into zones, with the same arithmetic as the
former iowave functions, so results stay bitwise identical, and adds the
zones given directly (ctrl.txt, repeatable; read into lexer, control.h):

  B 520 id method priority     method 1: relaxation (wave generation), 2: beach
  B 521 id edge s0 s1 width    edge 1: x-, 2: x+, 3: y-, 4: y+; along-edge range
                               [s0,s1] from the edge's start (s1 <= s0: whole edge)
  B 524 id source              source of the zone (repeatable; 1: the B 92 wave)
--------------------------------------------------------------------*/

enum class bc_method {relax, beach};

class bc_zone
{
public:
    bc_zone(int, bc_method, double, double, double, double, double, double);

    bool inside(double, double) const;   // in one of the two triangles of the quadrilateral
    double line_dist(double, double) const;  // distance to the reference line

    // the 4 numbers xs,xe,ys,ye of the axis-aligned box around the quadrilateral
    void box(double&, double&, double&, double&) const;

    int id;
    int priority = 0;            // where zones overlap, the highest priority sets the target
    bool user = false;           // given by B 520 (not translated from B 96 / 107 / 108)
    std::vector<int> sources;    // B 524 source ids; empty: all sources
    bc_method method;
    double xs,ys,xe,ye,d;   // reference line and half width
    double fac;             // beach: distance factor (2 for B 99 1)
    double P1[2],P2[2],P3[2],P4[2];

    static int intriangle(double,double,double,double,double,double,double,double);
};

class bc_zone_set
{
public:
    static bc_zone_set from_legacy(lexer*, ghostcell*);
    
    // the relaxation zone that sets the target at (x,y): highest priority, then
    // first given; nullptr outside all relaxation zones
    const bc_zone* relax_zone_at(double, double) const;
    
    bool has_sources() const;    // any zone with B 524 sources
    bool user_relax() const;     // any B 520 relaxation zone
    bool user_beach() const;     // any B 520 beach zone

    double relax_weight(double, double) const;   // generation zones (former rb1_ext)
    int relax_flag(double, double) const;        // inside a generation zone (former rb1_flag)
    double beach_weight(double, double) const;   // beach zones (former rb3_ext)
    double relax_dist(double, double) const;     // former distgen_calc
    double beach_dist(double, double) const;     // former distbeach_calc

    // no-refinement boxes for reefamr (xs,xe,ys,ye per box): the former B 96
    // ranges for old input (unchanged), plus the box of every B 520 zone
    void norefine_boxes(lexer*, std::vector<double>&) const;

    std::vector<bc_zone> relax, beach;
    
private:
    void read_input(lexer*, ghostcell*);
};

#endif
