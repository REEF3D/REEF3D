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

// The positions of the neighbours' particles within 3 cells of the subdomain boundary (the
// kernel of the free surface, pic_f_phi.cpp; the velocities to the grid, pic_f_flip.cpp): x first,
// then y with what came from x, then z with all, so those of the edge and corner neighbours arrive
// too. H: x y z of each, with vel also u v w.
void pic_f::halo(lexer *p, ghostcell *pgc, std::vector<double> &H, bool vel)
{
    const int band = 3;
    const int st = vel ? 15 : 3;   // with velocities: and the gradients of APIC
    H.clear();

    for(int dir=0; dir<3; ++dir)
    {
        if(dir==1 && p->j_dir==0)
        continue;

        const double *N = dir==0 ? p->XN : (dir==1 ? p->YN : p->ZN);
        const int nc = dir==0 ? p->knox : (dir==1 ? p->knoy : p->knoz);
        const int nblo = dir==0 ? p->nb1 : (dir==1 ? p->nb3 : p->nb5);
        const int nbhi = dir==0 ? p->nb4 : (dir==1 ? p->nb2 : p->nb6);

        std::vector<double> lo, hi, recv;

        auto pick = [&](const double *r)
        {
            const double c = r[dir==0 ? 0 : (dir==1 ? 1 : 2)];
            const int ci = cell(N,nc,c);

            if(ci<band && nblo!=-2)
            lo.insert(lo.end(),r,r+st);

            if(ci>=nc-band && nbhi!=-2)
            hi.insert(hi.end(),r,r+st);
        };

        for(const auto &q : P)
        {
        const double r[15] = {q.x,q.y,q.z,q.u,q.v,q.w,q.b[0],q.b[1],q.b[2],q.b[3],q.b[4],q.b[5],q.b[6],q.b[7],q.b[8]};
        pick(r);
        }

        const size_t nh = H.size();
        for(size_t m=0; m+st<=nh; m+=st)
        pick(&H[m]);

        sendrecv(pgc,nblo,nbhi,lo,hi,recv,8200+10*dir);

        H.insert(H.end(),recv.begin(),recv.end());
    }
}

// lo to the neighbour below (nblo), hi to the one above (nbhi); recv: what both sent here
void pic_f::sendrecv(ghostcell *pgc, int nblo, int nbhi, std::vector<double> &lo, std::vector<double> &hi,
                     std::vector<double> &recv, int tag)
{
    int nlo = int(lo.size()), nhi = int(hi.size());
    int rlo = 0, rhi = 0;
    MPI_Request req[8];
    int nreq = 0;

    // counts: down (to nblo) with tag, up (to nbhi) with tag+1
    if(nblo!=-2)
    {
    MPI_Isend(&nlo,1,MPI_INT,nblo,tag,pgc->mpi_comm,&req[nreq++]);
    MPI_Irecv(&rlo,1,MPI_INT,nblo,tag+1,pgc->mpi_comm,&req[nreq++]);
    }
    if(nbhi!=-2)
    {
    MPI_Isend(&nhi,1,MPI_INT,nbhi,tag+1,pgc->mpi_comm,&req[nreq++]);
    MPI_Irecv(&rhi,1,MPI_INT,nbhi,tag,pgc->mpi_comm,&req[nreq++]);
    }
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    std::vector<double> blo(size_t(rlo)+1), bhi(size_t(rhi)+1);
    lo.push_back(0.0);
    hi.push_back(0.0);
    nreq = 0;

    if(nblo!=-2)
    {
    MPI_Isend(lo.data(),nlo,MPI_DOUBLE,nblo,tag+2,pgc->mpi_comm,&req[nreq++]);
    MPI_Irecv(blo.data(),rlo,MPI_DOUBLE,nblo,tag+3,pgc->mpi_comm,&req[nreq++]);
    }
    if(nbhi!=-2)
    {
    MPI_Isend(hi.data(),nhi,MPI_DOUBLE,nbhi,tag+3,pgc->mpi_comm,&req[nreq++]);
    MPI_Irecv(bhi.data(),rhi,MPI_DOUBLE,nbhi,tag+2,pgc->mpi_comm,&req[nreq++]);
    }
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    recv.assign(blo.begin(),blo.begin()+rlo);
    recv.insert(recv.end(),bhi.begin(),bhi.begin()+rhi);
}
