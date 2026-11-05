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

#ifndef SANDSLIDE_H_
#define SANDSLIDE_H_

#include"sediment_curv.h"

class lexer;
class ghostcell;
class sediment_fdm;
class slice;
class sliceint;

using namespace std;

// sand slide transfer helpers (used inside the slide loops, with the increment macros IP, JP):
// SLIDE_AR: area ratio donor/receiver, a slid height a of the donor is a*SLIDE_AR on the receiver
//           (volume conserving on non-uniform grids)
// SLIDE_NB: the neighbour lies inside the global domain (no transfer into physical boundary ghost
//           cells, where it was lost) and, in 2D (j_dir 0), in the same row
#define SLIDE_AR(di,dj) ((p->DXN[IP]*p->DYN[JP])/(p->DXN[IP+(di)]*p->DYN[JP+(dj)]))
// curvilinear grid (G 50, the *_curv functions): SLIDE_ARC the area ratio of the geometry layer, SLIDE_D the
// distance between the centres of the cell and its neighbour (cart: the Cartesian expression), sediment_curv.h
#define SLIDE_ARC(di,dj) (sediment_curv::aratio(p,i,j,(di),(dj)))
#define SLIDE_D(di,dj,cart) (p->geo_curv==1 ? sediment_curv::dist(p,i,j,(di),(dj)) : (cart))
#define SLIDE_NB(di,dj) (i+p->origin_i+(di)>=0 && i+p->origin_i+(di)<p->gknox && j+p->origin_j+(dj)>=0 && j+p->origin_j+(dj)<p->gknoy && ((dj)==0 || (p->j_dir==1 && p->gknoy>1)))

// mesh refinement and multi-block grids (nhflow_amr, sflow_amr): only the cells with amr_donor > 0.5 give
// material (the cells of the grid; not the cells around it, the ring of a native grid, the covered cells)
#define SLIDE_DONOR if(amr_donor==nullptr || (*amr_donor)(i,j)>0.5)

class sandslide  
{
public:

    virtual ~sandslide(){}
	virtual void start(lexer*,ghostcell*,sediment_fdm*)=0;

    // mesh refinement and multi-block grids: one iteration of start() in steps, run by the driver over all
    // grids in turn, with the exchanges between the grids between the steps (sediment_amr.cpp)
    //   amr_begin: slide_fh = 0
    //   amr_sweep: fh = 0, the transfers of the donor cells into fh of the cell and its neighbours (also the
    //              ghost cells: the driver adds them to the cells they stand for); returns the number of
    //              over-steep cells
    //   amr_apply: slide_fh += fh, bedzh += fh in the sediment cells
    void amr_begin(lexer*,sediment_fdm*);
    virtual int amr_sweep(lexer*,ghostcell*,sediment_fdm*) { return 0; }
    void amr_apply(lexer*,sediment_fdm*);
    virtual slice* amr_fh() { return nullptr; }     // the transfers of the last sweep (nullptr: no slide)
    slice *amr_donor = nullptr;
};

#endif

