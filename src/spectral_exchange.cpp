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

#include"spectral_exchange.h"
#include"spectral_store.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<mpi.h>

// side s = 0..3: i- (nb1, gcslpara1), i+ (nb4, gcslpara4), j- (nb3, gcslpara3), j+ (nb2, gcslpara2)

static int nb_of(lexer *p, int s)
{
    if(s==0) return p->nb1;
    if(s==1) return p->nb4;
    if(s==2) return p->nb3;
    return p->nb2;
}

static int count_of(lexer *p, int s)
{
    if(s==0) return p->gcslpara1_count;
    if(s==1) return p->gcslpara4_count;
    if(s==2) return p->gcslpara3_count;
    return p->gcslpara2_count;
}

static int *cell_of(lexer *p, int s, int q)
{
    if(s==0) return p->gcslpara1[q];
    if(s==1) return p->gcslpara4[q];
    if(s==2) return p->gcslpara3[q];
    return p->gcslpara2[q];
}

// interior cell of layer n next to the boundary cell (i,j), and the halo cell it fills on the other side
static void interior(int s, int i, int j, int n, int &ii, int &jj)
{
    ii=i; jj=j;
    if(s==0) ii=i+n;
    if(s==1) ii=i-n;
    if(s==2) jj=j+n;
    if(s==3) jj=j-n;
}

static void halo(int s, int i, int j, int n, int &ii, int &jj)
{
    ii=i; jj=j;
    if(s==0) ii=i-n-1;
    if(s==1) ii=i+n+1;
    if(s==2) jj=j-n-1;
    if(s==3) jj=j+n+1;
}

spectral_exchange::spectral_exchange(lexer *p, int nbin_, int layers_) : nbin(nbin_), layers(std::min(std::max(layers_,1),p->margin)), sent(0)
{
    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0)
    {
    sendbuf[s].assign(size_t(count_of(p,s))*layers*nbin,0.0f);
    recvbuf[s].assign(size_t(count_of(p,s))*layers*nbin,0.0f);
    }
}

void spectral_exchange::pack(lexer *p, spectral_store &N, int s)
{
    size_t c=0;

    for(int q=0; q<count_of(p,s); ++q)
    {
    const int *ij = cell_of(p,s,q);

        for(int n=0; n<layers; ++n)
        {
        int ii,jj;
        interior(s,ij[0],ij[1],n,ii,jj);
        const float *f = N.spec(ii,jj);

            if(f!=nullptr && N.active(ii,jj))
            std::copy(f,f+nbin,&sendbuf[s][c]);
            else
            std::fill(&sendbuf[s][c],&sendbuf[s][c]+nbin,0.0f);

        c+=nbin;
        }
    }
}

void spectral_exchange::unpack(lexer *p, spectral_store &N, int s)
{
    size_t c=0;

    for(int q=0; q<count_of(p,s); ++q)
    {
    const int *ij = cell_of(p,s,q);

        for(int n=0; n<layers; ++n)
        {
        int ii,jj;
        halo(s,ij[0],ij[1],n,ii,jj);
        float *f = N.spec(ii,jj);

            if(f!=nullptr && N.active(ii,jj))
            std::copy(&recvbuf[s][c],&recvbuf[s][c]+nbin,f);

        c+=nbin;
        }
    }
}

void spectral_exchange::start(lexer *p, ghostcell *pgc, spectral_store &N)
{
    MPI_Request req[8];
    int nreq=0;

    // tag: the side of the sender; the receiver expects the opposite side
    static const int opposite[4] = {1,0,3,2};

    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0 && !recvbuf[s].empty())
    MPI_Irecv(recvbuf[s].data(),int(recvbuf[s].size()),MPI_FLOAT,nb_of(p,s),900+opposite[s],pgc->mpi_comm,&req[nreq++]);

    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0 && !sendbuf[s].empty())
    {
    pack(p,N,s);
    MPI_Isend(sendbuf[s].data(),int(sendbuf[s].size()),MPI_FLOAT,nb_of(p,s),900+s,pgc->mpi_comm,&req[nreq++]);
    sent += long(sendbuf[s].size()*sizeof(float));
    }

    if(nreq>0)
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0 && !recvbuf[s].empty())
    unpack(p,N,s);
}
