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

//  CFD in checkpoints (P 48, I 45): every field and slice of the level-0 grid (ckpt_registry.h), in
//  the order of their construction, so the arrays of all objects (fdm, momentum stages, turbulence,
//  level set, heat, concentration, ...) are covered.  The kinds of the arrays are kept with them: a
//  restart with other objects (other keys) is refused.  The work arrays of the solvers (local to the
//  rank) are restored on the same decomposition only; the scalars of the objects (raw arrays, the same
//  on all ranks) on any.

#include"cfd_ckpt.h"
#include"lexer.h"
#include"reefckpt.h"
#include"ckpt_registry.h"
#include<mpi.h>
#include<vector>
#include<string>
#include<iostream>

using namespace std;

namespace
{
vector<ckpt_registry::entry> grid_arrays(const void *owner)
{
    vector<ckpt_registry::entry> e;
    ckpt_registry::visit(owner,[&](const ckpt_registry::entry &x){ e.push_back(x); });
    return e;
}

double code(const ckpt_registry::entry &e)
{
    return (e.kind=='C' ? 0.0 : (e.kind=='S' ? 2.0 : 8.0)) + (e.isint ? 1.0 : 0.0) + (e.local ? 4.0 : 0.0);
}
}

void cfd_ckpt_write(lexer *p, reefckpt &ck, const string &tag)
{
    const vector<ckpt_registry::entry> e = grid_arrays(p);

    vector<double> kinds;
    for(const auto &x : e)
    kinds.push_back(code(x));
    ck.rank_blob(tag+"_kinds",kinds);

    for(size_t n=0; n<e.size(); ++n)
    {
        const string nm = tag+to_string(n);
        if(e[n].kind=='R')
        {
            ck.rank_blob(nm,vector<double>((const double*)e[n].V,(const double*)e[n].V+e[n].n));
            continue;
        }
        if(e[n].kind=='C')
        {
            if(e[n].isint) ck.put_Ci(nm,(const int*)e[n].V);
            else           ck.put_C(nm,(const double*)e[n].V);
        }
        else
        {
            if(e[n].isint) ck.put_Si(nm,(const int*)e[n].V);
            else           ck.put_S(nm,(const double*)e[n].V);
        }
    }
}

bool cfd_ckpt_read(lexer *p, reefckpt &ck, string &msg, const string &tag)
{
    const vector<ckpt_registry::entry> e = grid_arrays(p);

    // the kinds of the checkpoint (the own piece on the same decomposition, else piece 0)
    const vector<double> kinds = ck.get_rank_blob(tag+"_kinds",ck.same_decomposition() ? p->mpirank : 0);

    int n = (int)min(e.size(),kinds.size());
    int bad = 0;
    for(int q=0; q<n; ++q)
    if(code(e[q])!=kinds[q])
    bad = 1;

    int nall=0, badall=0;
    MPI_Allreduce(&n,&nall,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&bad,&badall,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(badall || kinds.empty())
    {
        msg = "the arrays of the checkpoint do not match the objects of this run (other keys?)";
        return false;
    }

    bool ok = true;
    for(int q=0; q<nall; ++q)
    {
        // the state of a rank (solver work arrays): on the same decomposition only
        if(e[q].local && !ck.same_decomposition())
        continue;
        const string nm = tag+to_string(q);
        if(e[q].kind=='R')
        {
            // (the scalars the same on all ranks: those of piece 0 on another decomposition)
            const vector<double> v = ck.get_rank_blob(nm,ck.same_decomposition() ? p->mpirank : 0);
            if((long long)v.size()==e[q].n)
            std::copy(v.begin(),v.end(),(double*)e[q].V);
            continue;
        }
        if(e[q].kind=='C')
        ok = (e[q].isint ? ck.get_Ci(nm,(int*)e[q].V) : ck.get_C(nm,(double*)e[q].V)) && ok;
        else
        ok = (e[q].isint ? ck.get_Si(nm,(int*)e[q].V) : ck.get_S(nm,(double*)e[q].V)) && ok;
    }

    if(p->mpirank==0 && (e.size()!=kinds.size()))
    cout<<"I 45: "<<kinds.size()<<" arrays in the checkpoint, "<<e.size()<<" in this run; "<<nall<<" restored"<<endl;

    return ok;
}
