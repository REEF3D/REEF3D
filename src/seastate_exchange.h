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

#ifndef SEASTATE_EXCHANGE_H_
#define SEASTATE_EXCHANGE_H_

#include"increment.h"
#include<vector>
#include<functional>

class lexer;
class ghostcell;
class seastate_store;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - halo exchange of the spectra between MPI ranks

Uses the 2D parallel boundary cell lists of the grid (gcslpara1-4, as
ghostcell::gcslparax) and the neighbour ranks nb1 (i-), nb4 (i+),
nb3 (j-), nb2 (j+), with its own buffers and non-blocking point-to-point
messages on pgc->mpi_comm. ghostcell itself is not changed.

  layers  number of halo layers (1 for first-order upwind, <= margin)

Cells without storage send zeros and ignore received values, so both
sides always exchange the same message sizes.
--------------------------------------------------------------------*/

// Phase 14: the second halo layer at the rank boundaries takes nv values per cell from the neighbour rank (its second
// interior layer): get(i,j,v) reads them from an interior cell, set(i,j,v) writes them into a halo cell. The grid flags
// of REEF3D mark only the first halo layer as part of the domain; the second-order fluxes (A 796 2) read two
void seastate_halo2(lexer*, ghostcell*, int nv, const std::function<void(int,int,double*)> &get,
                    const std::function<void(int,int,const double*)> &set);

// Phase 16: the first halo layer, nv values per cell, with the side s of the rank boundary (0 i-, 1 i+, 2 j-, 3 j+):
// get(s,i,j,v) reads them from the interior cell (i,j) next to side s, set(s,i,j,v) writes them into the halo cell
// (i,j) of side s
void seastate_halo1(lexer*, ghostcell*, int nv, int tag, const std::function<void(int,int,int,double*)> &get,
                    const std::function<void(int,int,int,const double*)> &set);

// Phase 16: nv floats on chosen cells of the first layer: send[s] the interior cells (i0,j0,i1,j1,...) next to side s,
// in the order of the boundary cells of the side, recv[s] the halo cells of side s that take the values of the
// neighbour rank in the same order (the neighbour's send list of the opposite side)
void seastate_halo_list(lexer*, ghostcell*, int nv, int tag, const vector<int> send[4], const vector<int> recv[4],
                        const std::function<void(int,int,int,float*)> &get, const std::function<void(int,int,int,const float*)> &set);

class seastate_exchange : public increment
{
public:
    seastate_exchange(lexer*, int nbin, int layers);

    void start(lexer*, ghostcell*, seastate_store&);

    long bytes_sent() const {return sent;}
    int halo_layers() const {return layers;}

private:
    void pack(lexer*, seastate_store&, int side);
    void unpack(lexer*, seastate_store&, int side);

    int nbin, layers;
    long sent;
    vector<float> sendbuf[4], recvbuf[4];
};

#endif
