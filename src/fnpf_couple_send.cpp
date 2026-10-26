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

#include"fnpf_couple_send.h"
#include"reefcouple.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<algorithm>

fnpf_couple_send::fnpf_couple_send(lexer *pp, fdm_fnpf *c, ghostcell *ppgc) : p(pp), pgc(ppgc)
{
    nz = p->knoz+1;
    
    i0g = p->origin_i;
    i1g = p->origin_i + p->knox;
    j0g = p->origin_j;
    j1g = p->origin_j + p->knoy;
    
    if(reefcouple::depth>8)
    {
    cout<<"fnpf_couple_send: depth > 8"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    L.resize(reefcouple::recv_ic.size());
    
    for(size_t n=0; n<L.size(); ++n)
    {
        L[n].ic  = reefcouple::recv_ic[n];
        L[n].app = reefcouple::recv_app[n];
        handshake(p,c,pgc,L[n]);
    }
    
    // level 0: the initial state
    send(p,c,pgc);
    
    if(p->mpirank==0)
    cout<<"FNPF coupling: sending to "<<L.size()<<" model(s), "<<nz<<" sigma nodes"<<endl;
}

fnpf_couple_send::~fnpf_couple_send()
{
}

void fnpf_couple_send::handshake(lexer *p, fdm_fnpf *c, ghostcell *pgc, link &lk)
{
    const int gx = p->gknox;
    const int gy = p->gknoy;
    const int m = increment::marge;   // offset of the 1D coordinate arrays (XP[IP])
    int i,j,k;
    
    // global cell centre coordinates
    std::vector<double> xs(gx,0.0), xc(gx,0.0), ys(gy,0.0), yc(gy,0.0);
    
    for(i=0; i<p->knox; ++i)
    {
    xs[p->origin_i+i] += p->XP[i+m];
    xc[p->origin_i+i] += 1.0;
    }
    
    for(j=0; j<p->knoy; ++j)
    {
    ys[p->origin_j+j] += p->YP[j+m];
    yc[p->origin_j+j] += 1.0;
    }
    
    MPI_Allreduce(MPI_IN_PLACE,xs.data(),gx,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,xc.data(),gx,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,ys.data(),gy,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,yc.data(),gy,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    
    // header to all ranks of the receiver
    const int root = p->mpirank==0 ? MPI_ROOT : MPI_PROC_NULL;
    
    int hi[8] = {reefcouple::version, gx, gy, nz, p->j_dir ? 1 : 0, p->mpi_size, reefcouple::depth, 0};
    MPI_Bcast(hi,8,MPI_INT,root,lk.ic);
    
    std::vector<double> hd(2+gx+gy+nz);
    hd[0] = p->wd;
    hd[1] = p->simtime;
    for(i=0; i<gx; ++i)
    hd[2+i] = xs[i]/std::max(xc[i],1.0);
    for(j=0; j<gy; ++j)
    hd[2+gx+j] = ys[j]/std::max(yc[j],1.0);
    for(k=0; k<nz; ++k)
    hd[2+gx+gy+k] = p->ZN[k+m];
    
    MPI_Bcast(hd.data(),int(hd.size()),MPI_DOUBLE,root,lk.ic);
    
    // index boxes: mine to the receiver, the receiver windows to me
    int nremote;
    MPI_Comm_remote_size(lk.ic,&nremote);
    
    int mine[4] = {i0g,i1g,j0g,j1g};
    std::vector<int> win(4*nremote);
    MPI_Allgather(mine,4,MPI_INT,win.data(),4,MPI_INT,lk.ic);
    
    lk.d.clear();
    for(int r=0; r<nremote; ++r)
    {
        dest dd;
        dd.rank = r;
        dd.i0 = std::max(i0g,win[4*r]);
        dd.i1 = std::min(i1g,win[4*r+1]);
        dd.j0 = std::max(j0g,win[4*r+2]);
        dd.j1 = std::min(j1g,win[4*r+3]);
        dd.ncol = (dd.i1-dd.i0)*(dd.j1-dd.j0);
        
        if(dd.i1>dd.i0 && dd.j1>dd.j0)
        lk.d.push_back(dd);
    }
    
    const int nd = int(lk.d.size());
    for(int s=0; s<reefcouple::depth; ++s)
    {
        lk.buf[s].resize(nd);
        lk.req[s].assign(nd,MPI_REQUEST_NULL);
        for(int q=0; q<nd; ++q)
        lk.buf[s][q].resize(1+size_t(lk.d[q].ncol)*(1+3*nz));
    }
    
    // still water depth: the bed, once
    std::vector<std::vector<double>> bb(nd);
    std::vector<MPI_Request> rq(nd);
    for(int q=0; q<nd; ++q)
    {
        const dest &dd = lk.d[q];
        bb[q].resize(dd.ncol);
        int n=0;
        for(int ig=dd.i0; ig<dd.i1; ++ig)
        for(int jg=dd.j0; jg<dd.j1; ++jg)
        bb[q][n++] = c->bed(ig-p->origin_i,jg-p->origin_j);
        
        MPI_Isend(bb[q].data(),dd.ncol,MPI_DOUBLE,dd.rank,reefcouple::tag_bed,lk.ic,&rq[q]);
    }
    MPI_Waitall(nd,rq.data(),MPI_STATUSES_IGNORE);
    
    int ntot=0, nlocal=nd;
    MPI_Reduce(&nlocal,&ntot,1,MPI_INT,MPI_SUM,0,pgc->mpi_comm);
    if(p->mpirank==0)
    cout<<"FNPF coupling: app "<<lk.app<<", "<<nremote<<" receiver ranks, "<<ntot<<" FNPF -> receiver messages per step"<<endl;
}

void fnpf_couple_send::pack(lexer *p, fdm_fnpf *c, const dest &dd, std::vector<double> &b)
{
    int i,j,k;
    const size_t nc = dd.ncol;
    
    b[0] = p->simtime;
    
    size_t n=0;
    for(int ig=dd.i0; ig<dd.i1; ++ig)
    for(int jg=dd.j0; jg<dd.j1; ++jg)
    {
        i = ig-p->origin_i;
        j = jg-p->origin_j;
        
        b[1+n] = c->eta(i,j);
        
        for(k=0; k<nz; ++k)
        {
        b[1+nc+n*nz+k]        = c->U[FIJK];
        b[1+nc+(nc+n)*nz+k]   = c->V[FIJK];
        b[1+nc+(2*nc+n)*nz+k] = c->W[FIJK];
        }
        ++n;
    }
}

void fnpf_couple_send::send(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const int s = level%reefcouple::depth;
    
    for(link &lk : L)
    if(lk.open)
    {
        // at most depth levels ahead of the receiver
        MPI_Waitall(int(lk.req[s].size()),lk.req[s].data(),MPI_STATUSES_IGNORE);
        
        for(size_t q=0; q<lk.d.size(); ++q)
        {
        pack(p,c,lk.d[q],lk.buf[s][q]);
        MPI_Issend(lk.buf[s][q].data(),int(lk.buf[s][q].size()),MPI_DOUBLE,lk.d[q].rank,reefcouple::tag_data,lk.ic,&lk.req[s][q]);
        }
    }
    
    ++level;
    
    check_stop(p,pgc);
}

void fnpf_couple_send::check_stop(lexer *p, ghostcell *pgc)
{
    const int nl = int(L.size());
    std::vector<int> stop(nl,0);
    
    if(p->mpirank==0)
    for(int n=0; n<nl; ++n)
    if(L[n].open)
    {
        int flag=0;
        MPI_Iprobe(0,reefcouple::tag_stop,L[n].ic,&flag,MPI_STATUS_IGNORE);
        if(flag)
        {
        int v;
        MPI_Recv(&v,1,MPI_INT,0,reefcouple::tag_stop,L[n].ic,MPI_STATUS_IGNORE);
        stop[n]=1;
        }
    }
    
    MPI_Bcast(stop.data(),nl,MPI_INT,0,pgc->mpi_comm);
    
    for(int n=0; n<nl; ++n)
    if(stop[n]==1)
    {
        close(L[n]);
        if(p->mpirank==0)
        cout<<"FNPF coupling: app "<<L[n].app<<" has stopped"<<endl;
    }
}

void fnpf_couple_send::close(link &lk)
{
    if(!lk.open)
    return;
    
    // end markers after the last level, then all sends complete (the receiver drains)
    lk.endval = p->simtime;
    lk.endreq.assign(lk.d.size(),MPI_REQUEST_NULL);
    for(size_t q=0; q<lk.d.size(); ++q)
    MPI_Isend(&lk.endval,1,MPI_DOUBLE,lk.d[q].rank,reefcouple::tag_end,lk.ic,&lk.endreq[q]);
    
    for(int s=0; s<reefcouple::depth; ++s)
    MPI_Waitall(int(lk.req[s].size()),lk.req[s].data(),MPI_STATUSES_IGNORE);
    
    MPI_Waitall(int(lk.endreq.size()),lk.endreq.data(),MPI_STATUSES_IGNORE);
    
    lk.open = false;
}

bool fnpf_couple_send::done() const
{
    for(const link &lk : L)
    if(lk.open)
    return false;
    
    return !L.empty();
}

void fnpf_couple_send::finish()
{
    if(finished)
    return;
    
    for(link &lk : L)
    if(lk.open)
    {
        close(lk);
        
        // the receiver sends its stop when it ends
        if(p->mpirank==0)
        {
        int v;
        MPI_Recv(&v,1,MPI_INT,0,reefcouple::tag_stop,lk.ic,MPI_STATUS_IGNORE);
        }
    }
    
    finished = true;
}
