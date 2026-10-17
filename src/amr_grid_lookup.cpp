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

#include"amr_grid_lookup.h"
#include<mpi.h>
#include<vector>

void amr_point_owners(const int *lev, int n, char *mine)
{
    if(n<=0)
    return;

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // the key of a point: the level first, the lower rank before the higher one
    std::vector<int> key(n), top(n);
    for(int m=0; m<n; ++m)
    key[m] = (lev[m]<0) ? -1 : lev[m]*size + (size-1-rank);
    MPI_Allreduce(key.data(),top.data(),n,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    for(int m=0; m<n; ++m)
    {
        if(top[m]<size)     // level 0 at most: every rank that holds the point
        mine[m] = (lev[m]==0) ? 1 : 0;
        else
        mine[m] = (key[m]==top[m]) ? 1 : 0;
    }
}

void amr_points::find(const amr_column_lookup *amr, lexer *p, int n, const double *x, const double *y, const int *flag)
{
    lev.assign(n,-1);
    q.assign(n,p);
    f.assign(n,nullptr);
    mine.assign(n,0);
    for(int m=0; m<n; ++m)
    {
        lexer *qq = p;
        void *ff = nullptr;
        const int L = amr->finest_column(x[m],y[m],qq,ff);
        if(L>0)
        {
            lev[m] = L;
            q[m] = qq;
            f[m] = ff;
        }
        else
        lev[m] = (flag[m]>0) ? 0 : -1;
    }
    amr_point_owners(lev.data(),n,mine.data());
}

void amr_cfd_mine(const amr_grid_lookup *amr, int n, const double *x, const double *y, const double *z, const int *flag, std::vector<char> &mine)
{
    mine.assign(n,0);
    if(amr==nullptr || !amr->placed())
    {
        for(int m=0; m<n; ++m)
        mine[m] = flag[m]>0 ? 1 : 0;
        return;
    }
    std::vector<int> lev(n,-1);
    for(int m=0; m<n; ++m)
    {
        const int l = amr->finest_level(x[m],y[m],z[m]);
        lev[m] = (l>0) ? l : (flag[m]>0 ? 0 : -1);
    }
    amr_point_owners(lev.data(),n,mine.data());
}
