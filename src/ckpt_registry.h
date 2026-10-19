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

#ifndef CKPT_REGISTRY_H_
#define CKPT_REGISTRY_H_

//  The live fields and slices of a run (field_base, slice_base), in the order of their construction,
//  for the checkpoints of CFD (P 48, I 45; driver_ckpt.cpp): every array a grid keeps over the steps
//  is in a field or a slice of the objects of that grid, so the checkpoint holds them all without a
//  list per class.  owner: the lexer of the grid (level 0, or a patch of the mesh refinement).
//  Other arrays of a grid join with add/remove: the work arrays of the iterative solvers (whose ghost
//  cells carry over from one solve to the next) as raw arrays local to the rank, the scalars of an
//  object that are the same on all ranks (a reference volume, a flag) as raw arrays of the grid.

#include<map>
#include<mutex>
#include<functional>

class ckpt_registry
{
public:
    struct entry
    {
        const void *owner;
        char kind;          // 'C' field (imax jmax kmax), 'S' slice (imax jmax), 'R' raw array
        bool isint;
        void *V;
        bool local;         // the state of a rank (solver work arrays): restored on the same decomposition only
        long long n = 0;    // 'R' raw array: its length (local: a matrix, a vector of the solver; else
                            // scalars the same on all ranks, restored on any decomposition)
    };

    static unsigned long long add(const void *owner, char kind, bool isint, void *V, bool local=false);
    static unsigned long long add_raw(const void *owner, double *V, long long n, bool local=true);
    static void remove(unsigned long long id);

    // the entries of one owner in the order of construction
    static void visit(const void *owner, const std::function<void(const entry&)>&);

private:
    static std::map<unsigned long long,entry>& table();
    static std::mutex& lock();
};

#endif
