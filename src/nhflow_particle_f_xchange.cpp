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

#include"nhflow_particle_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<cstring>

// Particles that left the subdomain are sent to the face neighbour (x first, then y).
// Particles crossing a corner arrive at the x-neighbour and are forwarded in the next pass.
// direction q: 0 -x (nb1), 1 +y (nb2), 2 -y (nb3), 3 +x (nb4)

void nhflow_particle_f::xchange(lexer *p, ghostcell *pgc)
{
    const int nb[4] = {p->nb1, p->nb2, p->nb3, p->nb4};
    const int opp[4] = {3,2,1,0};
    const int tagbase = 7700;

    for(int pass=0; pass<6; ++pass)
    {
        std::vector<double> sendbuf[4], recvbuf[4];
        int sendnum[4]={0,0,0,0}, recvnum[4]={0,0,0,0};

        size_t q=0;
        while(q<P.size())
        {
            const nhflow_particle_data &a = P[q];
            int dir=-1;

            if(a.x<xloc_s && nb[0]!=-2)
            dir=0;
            else if(a.x>=xloc_e && nb[3]!=-2)
            dir=3;
            else if(p->j_dir==1 && a.y<yloc_s && nb[2]!=-2)
            dir=2;
            else if(p->j_dir==1 && a.y>=yloc_e && nb[1]!=-2)
            dir=1;

            if(dir>=0)
            {
                const double *v = reinterpret_cast<const double*>(&a);
                sendbuf[dir].insert(sendbuf[dir].end(),v,v+NF);
                ++sendnum[dir];

                P[q] = P.back();
                P.pop_back();
                continue;
            }

            ++q;
        }

        int moved = sendnum[0]+sendnum[1]+sendnum[2]+sendnum[3];
        moved = pgc->globalisum(moved);

        if(moved==0)
        break;

        MPI_Request req[8];
        int nreq=0;

        for(int dir=0; dir<4; ++dir)
        if(nb[dir]!=-2)
        {
        MPI_Isend(&sendnum[dir],1,MPI_INT,nb[dir],tagbase+dir,pgc->mpi_comm,&req[nreq++]);
        MPI_Irecv(&recvnum[dir],1,MPI_INT,nb[dir],tagbase+opp[dir],pgc->mpi_comm,&req[nreq++]);
        }
        MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

        nreq=0;
        for(int dir=0; dir<4; ++dir)
        if(nb[dir]!=-2)
        {
        recvbuf[dir].resize(size_t(recvnum[dir])*NF+1);
        sendbuf[dir].resize(size_t(sendnum[dir])*NF+1);
        MPI_Isend(sendbuf[dir].data(),sendnum[dir]*NF,MPI_DOUBLE,nb[dir],tagbase+10+dir,pgc->mpi_comm,&req[nreq++]);
        MPI_Irecv(recvbuf[dir].data(),recvnum[dir]*NF,MPI_DOUBLE,nb[dir],tagbase+10+opp[dir],pgc->mpi_comm,&req[nreq++]);
        }
        MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

        for(int dir=0; dir<4; ++dir)
        for(int n2=0; n2<recvnum[dir]; ++n2)
        {
            nhflow_particle_data a;
            std::memcpy(&a,&recvbuf[dir][size_t(n2)*NF],sizeof(nhflow_particle_data));
            P.push_back(a);
        }
    }
}
