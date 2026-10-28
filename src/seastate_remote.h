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


#ifndef SEASTATE_REMOTE_H_
#define SEASTATE_REMOTE_H_

#include"seastate_coupling.h"
#include<mpi.h>

class lexer;
class ghostcell;

/*--------------------------------------------------------------------
REEF3D::SEASTATE remote waves (A 750 2): SEASTATE as its own model next
to the host (SFLOW A 10 2 or NHFLOW A 10 5) in one MPI job
(reef3d_mpmd), each with its own working directory, grid, bathymetry,
decomposition and number of ranks:

  mpirun -np 4 -wdir surge REEF3D : -np 2 -wdir waves REEF3D

  surge/ctrl.txt  A 10 2 (or 5), A 750 2, A 706, A 751, A 752, A 753
  waves/ctrl.txt  A 10 7, A 750 2, the same A 706, A 751, A 753, and the
                  wave physics, boundary (A 711) and wind (A 730 2,
                  W 200 ...) of a stand-alone run; the same model frame
                  and start date-time (W 203 / A 780) as the host

The waves step A 706 s at a time between the exchanges and do not stop
on their own: the host ends them at its end (ghostcell::final). The
wave grid is typically coarser than the host grid (the wave field is
smooth on the scale of the surge grid), so that the wave model costs
a fraction of the in-process coupling (A 750 1), and the two models run
at the same time. The exchange and its interpolation are in
seastate_coupling (remote_ini, remote_step) and seastate_remote_wave.

roles: collective on MPI_COMM_WORLD after lexer_read (driver), finds
the host and the wave app and creates their intercommunicator; the
rebalance stop (G 46) is switched off in both.
--------------------------------------------------------------------*/

class seastate_remote
{
public:
    static void roles(lexer*, ghostcell*);

    // called by ghostcell::final (normal end): the host stops the waves
    static void finish();

    static bool host() {return role==1;}
    static bool waves() {return role==2;}

    static MPI_Comm ic;                     // intercommunicator host <-> waves
    static int role;                        // 0 none, 1 host, 2 waves
    static seastate_coupling *link;         // host: the coupling object (stop message)
    static lexer *plex;

    static const int tag_ic = 4900;
};

// the wave side: SEASTATE on its own grid, driven by the host's exchanges
class seastate_remote_wave : public seastate_coupling
{
public:
    seastate_remote_wave(lexer*);
    virtual ~seastate_remote_wave() = default;

    // set-up, then wave steps until the host stops; ends with ghostcell::final
    void run(lexer*, ghostcell*);

private:
    // the host fields at the wave cells: the area mean over the host cells inside a wave cell (eta, U, V over the wet ones,
    // WL over all), bilinear where the host grid is coarser; cells outside the host grid keep the own still water
    void map_host(lexer*);
    bool ic_null() const;

    slice4 etaW, WLW, UW, VW, wetd;
    std::vector<double> hwx, hwy;           // widths of the host cells (area weights of the box mean)
};

#endif
