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
#include"ghostcell.h"
#include<mpi.h>
#include<cstring>

// Particles outside the subdomain go to the face neighbour on that side: x first, then y,
// then z. A particle that crossed an edge or a corner reaches the x-neighbour first and is
// passed on in the next pass; the passes stop when no particle moved on any rank.
// direction d: 0 -x (nb1), 1 +x (nb4), 2 -y (nb3), 3 +y (nb2), 4 -z (nb5), 5 +z (nb6);
// a rank receives from direction d what its neighbour there sent towards it (opp[d]).

void pic_f::xchange(lexer *p, ghostcell *pgc)
{
    const int nb[6] = {p->nb1, p->nb4, p->nb3, p->nb2, p->nb5, p->nb6};
    const int opp[6] = {1,0,3,2,5,4};
    const int tagbase = 8100;

    sent = 0;

    for(int pass=0; pass<8; ++pass)
    {
        std::vector<double> sendbuf[6], recvbuf[6];
        int sendnum[6]={0,0,0,0,0,0}, recvnum[6]={0,0,0,0,0,0};

        size_t n=0;
        while(n<P.size())
        {
            const pic_data &q = P[n];
            int d=-1;

            if(q.x<xloc_s && nb[0]!=-2)
            d=0;
            else if(q.x>=xloc_e && nb[1]!=-2)
            d=1;
            else if(p->j_dir==1 && q.y<yloc_s && nb[2]!=-2)
            d=2;
            else if(p->j_dir==1 && q.y>=yloc_e && nb[3]!=-2)
            d=3;
            else if(q.z<zloc_s && nb[4]!=-2)
            d=4;
            else if(q.z>=zloc_e && nb[5]!=-2)
            d=5;

            if(d>=0)
            {
                const double *v = reinterpret_cast<const double*>(&q);
                sendbuf[d].insert(sendbuf[d].end(),v,v+NF);
                ++sendnum[d];

                P[n] = P.back();
                P.pop_back();
                continue;
            }

            ++n;
        }

        int moved = sendnum[0]+sendnum[1]+sendnum[2]+sendnum[3]+sendnum[4]+sendnum[5];
        sent += moved;
        moved = pgc->globalisum(moved);

        if(moved==0)
        break;

        MPI_Request req[12];
        int nreq=0;

        for(int d=0; d<6; ++d)
        if(nb[d]!=-2)
        {
        MPI_Isend(&sendnum[d],1,MPI_INT,nb[d],tagbase+d,pgc->mpi_comm,&req[nreq++]);
        MPI_Irecv(&recvnum[d],1,MPI_INT,nb[d],tagbase+opp[d],pgc->mpi_comm,&req[nreq++]);
        }
        MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

        nreq=0;
        for(int d=0; d<6; ++d)
        if(nb[d]!=-2)
        {
        recvbuf[d].resize(size_t(recvnum[d])*NF+1);
        sendbuf[d].resize(size_t(sendnum[d])*NF+1);
        MPI_Isend(sendbuf[d].data(),sendnum[d]*NF,MPI_DOUBLE,nb[d],tagbase+10+d,pgc->mpi_comm,&req[nreq++]);
        MPI_Irecv(recvbuf[d].data(),recvnum[d]*NF,MPI_DOUBLE,nb[d],tagbase+10+opp[d],pgc->mpi_comm,&req[nreq++]);
        }
        MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

        for(int d=0; d<6; ++d)
        for(int m=0; m<recvnum[d]; ++m)
        {
            pic_data q;
            std::memcpy(&q,&recvbuf[d][size_t(m)*NF],sizeof(pic_data));
            P.push_back(q);
        }
    }

    sent = (long long)pgc->globalsum(double(sent));
    sent_total += sent;
}
