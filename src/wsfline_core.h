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

#ifndef WSFLINE_CORE_H_
#define WSFLINE_CORE_H_

#include<fstream>
#include<functional>

class lexer;
class ghostcell;

using namespace std;

// shared part of the free-surface line output of all solvers (P 52 lines along x, P 56 lines along y):
// local buffers, gather on rank 0, sort along the line, removal of double entries and the data rows.
// The solver adapters (print_wsfline_x/_y, nhflow_, fnpf_, sflow_print_wsfline/_y) locate the lines,
// fill loc / wsf / flag and write their file header.
class wsfline_core
{
public:
    // nline lines with nloc local points each (knox for lines along x, knoy for lines along y)
    wsfline_core(lexer*, ghostcell*, int nline, int nloc);
    ~wsfline_core();

    // loc = 1e20, wsf = -1e20 for all local points
    void reset();

    // gather the lines on rank 0 and write one row per point along the lines:
    // coordinate and value with setprecision(prec), coord(loc) is the coordinate written,
    // theory(loc) (if given) the wave-theory column; fill is written twice for a line without
    // a point in that row
    void write_rows(lexer*, ghostcell*, ofstream&, int prec, const function<double(double)> &coord,
                    const function<double(double)> &theory, const char *fill);

    double **loc, **wsf;
    int **flag;
    const int nline;
    int maxn, sumn;

private:
    void sort(double*, double*, int*, int, int);
    void remove_multientry(lexer*, double*, double*, int*, int&);

    double **loc_all, **wsf_all;
    int **flag_all, *rowflag, *wsfpoints;
};

#endif
