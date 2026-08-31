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
#define SLIDE_NB(di,dj) (i+p->origin_i+(di)>=0 && i+p->origin_i+(di)<p->gknox && j+p->origin_j+(dj)>=0 && j+p->origin_j+(dj)<p->gknoy && ((dj)==0 || (p->j_dir==1 && p->gknoy>1)))

class sandslide  
{
public:
    virtual ~sandslide() = default;

	virtual void start(lexer*,ghostcell*,sediment_fdm*)=0;
};

#endif

