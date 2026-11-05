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

#ifndef RIVER_INI_H_
#define RIVER_INI_H_

#include"increment.h"
#include<vector>

class lexer;
class ghostcell;
class slice;

using namespace std;

// River free-surface initialisation along the i-direction of the grid (F 65), shared by SFLOW and
// NHFLOW, for rivers whose water level follows the slope of the bed. The cross-sections are the grid
// lines i = const: the curvilinear corridor grid (G 50) has i along the river from the inflow; on a
// Cartesian grid the river runs along x.
//
//   F 65 1  water level linear in the distance along the river, from F 61 at the inflow column to
//           F 62 at the outflow column
//   F 65 2  steady 1D flow of the discharge W 10: standard-step backwater from F 62 at the outflow
//           column upstream, energy equation with the bed friction of B 50 and strip conveyance of the
//           cells of each cross-section (SFLOW: Manning n = ks^(1/6)/20 as sflow_rough_manning; NHFLOW:
//           the log law at the centre of the bottom layer as nhflow_bcmom A 519 1, depth-uniform
//           velocity); a column without a subcritical solution gets the critical level. The velocities of
//           the wet cells: strip velocity of the column friction slope along the i-direction of the cell,
//           so the initial discharge of every cross-section is W 10 (NHFLOW: on all layers).
//
// The level is uniform in a cross-section; cells above it, and side pools below it that are not
// connected to the lowest cell of the cross-section (neighbours in j), start dry. All ranks gather the
// global cross-sections (gknox x gknoy cells) and solve the same 1D problem.

class river_ini : public increment
{
public:
    river_ini();

    // solve for the levels from the bed levels of the cells; 0: F 65 off (or not on this grid), 1: done;
    // verbose 1 prints the levels at the inflow and outflow
    // law 0: Manning n = ks^(1/6)/20 (SFLOW); 1: log law at the centre of the bottom layer (NHFLOW A 519 1)
    int start(lexer*, ghostcell*, slice &bed, int verbose, int law);

    // eta of cell (ii,jj) (local indices) relative to p->wd: the column level for the connected wet
    // cells, the bed (depth 0) for the others
    double eta(lexer*, int ii, int jj, double bed) const;

    // F 65 2: Cartesian velocity of cell (ii,jj); 0 for dry cells
    void velocity(lexer*, int ii, int jj, double bed, double &u, double &v) const;

private:
    struct xsec
    {
        std::vector<double> zb, w, ds;     // bed, width across, length along of the cells of the column
        std::vector<int> J;                // their j index
        int nmin=0;                        // the lowest cell
    };

    std::vector<xsec> xs;
    std::vector<double> wl;
    std::vector<char> inrun;
    int NI=0, NJ=0;
    double Q=0.0, kst=0.0, g=9.81, ksl=0.05, dz0=1.0;
    int law=0;

    double strip(double d) const;

    static void xs_run(const xsec&, double h, int &n0, int &n1);
    static double xs_area(const xsec&, double h);
    double xs_conv(const xsec&, double h, double kst, double &dsw) const;
    static double xs_zmin(const xsec&);
    static double xs_dsmean(const xsec&);
    static double xs_energy(const xsec&, double h, double Q, double g);
    static double xs_critical(const xsec&, double Q, double g, double hmax);
};

#endif
