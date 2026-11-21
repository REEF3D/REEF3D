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

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<mpi.h>
#include<iostream>

// Every fluid cell (flag4 > 0) is split into F 83 sub-cells per direction (x, z; and y in 3D).
// A particle is put at the centre of every sub-cell that lies in the water (phi > 0, the level
// set interpolated to the sub-cell) and not in a solid (pic_f::blocked); for passive particles
// (F 81 10) up to a quarter of the sub-cell off the centre, from a hash of the particle id, so the
// seeding does not depend on the decomposition.

void pic_f::seed(lexer *p, fdm *a, ghostcell *pgc)
{
    P.clear();

    const int ny = p->j_dir==1 ? npd : 1;

    LOOP
    {
        // the global cell number and the first id of its particles
        const long long gi = (long long)(i + p->origin_i);
        const long long gj = (long long)(j + p->origin_j);
        const long long gk = (long long)(k + p->origin_k);
        const long long gcell = (gk*(long long)p->gknoy + gj)*(long long)p->gknox + gi;

        for(int qi=0; qi<npd; ++qi)
        for(int qj=0; qj<ny; ++qj)
        for(int qk=0; qk<npd; ++qk)
        {
            const int s = (qi*ny + qj)*npd + qk;
            const uint64_t id = uint64_t(gcell)*uint64_t(npc) + uint64_t(s);

            // F 81 10 (passive): up to a quarter of the sub-cell off its centre; F 81 1, 2: at the
            // centre (a random offset makes the free surface rough, and in still water the
            // roughness starts waves of a cell's length)
            const double jit = surface ? 0.0 : 0.5;
            const double ox = jit*(hash_uniform(id,0)-0.5);
            const double oy = jit*(hash_uniform(id,1)-0.5);
            const double oz = jit*(hash_uniform(id,2)-0.5);

            pic_data q;
            q.x = p->XN[IP] + (double(qi)+0.5+ox)/double(npd)*p->DXN[IP];
            q.y = p->j_dir==1 ? p->YN[JP] + (double(qj)+0.5+oy)/double(npd)*p->DYN[JP] : p->YP[JP];
            q.z = p->ZN[KP] + (double(qk)+0.5+oz)/double(npd)*p->DZN[KP];
            q.u = q.v = q.w = 0.0;
            q.id = double(id);

            if(ipol(p,a->phi,4,q.x,q.y,q.z)<=0.0)
            continue;

            if(blocked(p,a,q.x,q.y,q.z))
            continue;

            P.push_back(q);
        }
    }

    long long mine = (long long)P.size();
    MPI_Allreduce(&mine,&seeded,1,MPI_LONG_LONG,MPI_SUM,pgc->mpi_comm);

    const double maxid = double(p->gknox)*double(p->gknoy)*double(p->gknoz)*double(npc);

    if(p->mpirank==0)
    {
    cout<<"PIC seeding: "<<seeded<<" particles"<<endl;

    if(maxid>2147483647.0)
    cout<<"PIC: the particle ids exceed 2^31, the id in the output (int32) is not unique"<<endl;
    }
}

double pic_f::hash_uniform(uint64_t id, uint64_t c)
{
    uint64_t z = id*0x9E3779B97F4A7C15ULL + (c+1)*0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z>>30))*0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z>>27))*0x94D049BB133111EBULL;
    z = z ^ (z>>31);

    return double(z>>11)*(1.0/9007199254740992.0);
}
