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

#include"seastate_exchange.h"
#include"seastate_store.h"
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

seastate_exchange::seastate_exchange(lexer *p, int nbin_, int layers_) : nbin(nbin_), layers(std::min(std::max(layers_,1),p->margin)), sent(0)
{
    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0)
    {
    sendbuf[s].assign(size_t(count_of(p,s))*layers*nbin,0.0f);
    recvbuf[s].assign(size_t(count_of(p,s))*layers*nbin,0.0f);
    }
}

void seastate_exchange::pack(lexer *p, seastate_store &N, int s)
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

void seastate_exchange::unpack(lexer *p, seastate_store &N, int s)
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

void seastate_halo2(lexer *p, ghostcell *pgc, int nv, const std::function<void(int,int,double*)> &get,
                    const std::function<void(int,int,const double*)> &set)
{
    static const int opposite[4] = {1,0,3,2};
    vector<double> sb[4], rb[4];
    MPI_Request req[8];
    int nreq=0;

    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0 && count_of(p,s)>0)
    {
    const int n = count_of(p,s);
    sb[s].assign(size_t(n)*nv,0.0);
    rb[s].assign(size_t(n)*nv,0.0);
        for(int q=0; q<n; ++q)
        {
        const int *ij = cell_of(p,s,q);
        int ii,jj;
        interior(s,ij[0],ij[1],1,ii,jj);
        get(ii,jj,&sb[s][size_t(q)*nv]);
        }
    MPI_Irecv(rb[s].data(),int(rb[s].size()),MPI_DOUBLE,nb_of(p,s),940+opposite[s],pgc->mpi_comm,&req[nreq++]);
    MPI_Isend(sb[s].data(),int(sb[s].size()),MPI_DOUBLE,nb_of(p,s),940+s,pgc->mpi_comm,&req[nreq++]);
    }

    if(nreq>0)
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    for(int s=0; s<4; ++s)
    if(!rb[s].empty())
        for(int q=0; q<count_of(p,s); ++q)
        {
        const int *ij = cell_of(p,s,q);
        int ii,jj;
        halo(s,ij[0],ij[1],1,ii,jj);
        set(ii,jj,&rb[s][size_t(q)*nv]);
        }
}

void seastate_halo1(lexer *p, ghostcell *pgc, int nv, int tag, const std::function<void(int,int,int,double*)> &get,
                    const std::function<void(int,int,int,const double*)> &set)
{
    static const int opposite[4] = {1,0,3,2};
    vector<double> sb[4], rb[4];
    MPI_Request req[8];
    int nreq=0;

    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0 && count_of(p,s)>0)
    {
    const int n = count_of(p,s);
    sb[s].assign(size_t(n)*nv,0.0);
    rb[s].assign(size_t(n)*nv,0.0);
        for(int q=0; q<n; ++q)
        {
        const int *ij = cell_of(p,s,q);
        get(s,ij[0],ij[1],&sb[s][size_t(q)*nv]);
        }
    MPI_Irecv(rb[s].data(),int(rb[s].size()),MPI_DOUBLE,nb_of(p,s),tag+opposite[s],pgc->mpi_comm,&req[nreq++]);
    MPI_Isend(sb[s].data(),int(sb[s].size()),MPI_DOUBLE,nb_of(p,s),tag+s,pgc->mpi_comm,&req[nreq++]);
    }

    if(nreq>0)
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    for(int s=0; s<4; ++s)
    if(!rb[s].empty())
        for(int q=0; q<count_of(p,s); ++q)
        {
        const int *ij = cell_of(p,s,q);
        int ii,jj;
        halo(s,ij[0],ij[1],0,ii,jj);
        set(s,ii,jj,&rb[s][size_t(q)*nv]);
        }
}

void seastate_halo_list(lexer *p, ghostcell *pgc, int nv, int tag, const vector<int> send[4], const vector<int> recv[4],
                        const std::function<void(int,int,int,float*)> &get, const std::function<void(int,int,int,const float*)> &set)
{
    static const int opposite[4] = {1,0,3,2};
    vector<float> sb[4], rb[4];
    MPI_Request req[8];
    int nreq=0;

    for(int s=0; s<4; ++s)
    if(nb_of(p,s)>=0)
    {
        if(!recv[s].empty())
        {
        rb[s].assign(recv[s].size()/2*size_t(nv),0.0f);
        MPI_Irecv(rb[s].data(),int(rb[s].size()),MPI_FLOAT,nb_of(p,s),tag+opposite[s],pgc->mpi_comm,&req[nreq++]);
        }

        if(!send[s].empty())
        {
        const size_t n = send[s].size()/2;
        sb[s].assign(n*size_t(nv),0.0f);
        for(size_t q=0; q<n; ++q)
        get(s,send[s][2*q],send[s][2*q+1],&sb[s][q*nv]);
        MPI_Isend(sb[s].data(),int(sb[s].size()),MPI_FLOAT,nb_of(p,s),tag+s,pgc->mpi_comm,&req[nreq++]);
        }
    }

    if(nreq>0)
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    for(int s=0; s<4; ++s)
    for(size_t q=0; q<rb[s].size()/size_t(nv); ++q)
    set(s,recv[s][2*q],recv[s][2*q+1],&rb[s][q*nv]);
}

void seastate_exchange::start(lexer *p, ghostcell *pgc, seastate_store &N)
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
