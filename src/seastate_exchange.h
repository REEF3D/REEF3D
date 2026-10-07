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

class seastate_exchange : public increment
{
public:
    seastate_exchange(lexer*, int nbin, int layers);

    void start(lexer*, ghostcell*, seastate_store&);

    long bytes_sent() const {return sent;}

private:
    void pack(lexer*, seastate_store&, int side);
    void unpack(lexer*, seastate_store&, int side);

    int nbin, layers;
    long sent;
    vector<float> sendbuf[4], recvbuf[4];
};

#endif
