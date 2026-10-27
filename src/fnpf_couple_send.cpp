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
#include"6DOF_external.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice.h"
#include<iostream>
#include<algorithm>
#include<cmath>

// the FNPF bodies follow the 6DOF bodies of a receiver (I 243): sixdof_external without own dynamics,
// the pose and velocity extrapolated from the receiver's samples to the stage times of FNPF
class reefcouple_bodies final : public sixdof_external
{
public:
    reefcouple_bodies(const fnpf_couple_send *s, int nb) : S(s), nb_(nb), c0(nb,Eigen::Vector3d::Zero()), R0(nb,Eigen::Matrix3d::Identity()) {}
    
    bool drives(int n) const override {return n<nb_;}
    
    void connect(int n, double, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Matrix3d&) override
    {
        if(n<nb_)
        {
        c0[n] = c;
        R0[n] = R;
        }
    }
    
    bool predictor() const override {return false;}
    void fluid_loads(double, const std::vector<vec6>&, const std::vector<mat6>&, bool) override {}
    void advance(double, double) override {}
    
    void predict(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) override
    {
        state(n,t,c,R,v,w);
    }
    
    void state(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) override
    {
        double s[18];
        
        if(S->body_state(n,t,s))
        {
            for(int q=0; q<3; ++q)
            {
            c(q) = s[q];
            v(q) = s[3+q];
            w(q) = s[6+q];
            }
            for(int r=0; r<3; ++r)
            for(int q=0; q<3; ++q)
            R(r,q) = s[9+3*r+q];
        }
        else
        {
        c = c0[n];
        R = R0[n];
        v.setZero();
        w.setZero();
        }
    }
    
private:
    const fnpf_couple_send *S;
    int nb_;
    std::vector<Eigen::Vector3d> c0;
    std::vector<Eigen::Matrix3d> R0;
};

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
    
    // global cell centre coordinates
    const int gx = p->gknox;
    const int gy = p->gknoy;
    const int m = increment::marge;   // offset of the 1D coordinate arrays (XP[IP])
    std::vector<double> xs(gx,0.0), xc(gx,0.0), ys(gy,0.0), yc(gy,0.0);
    
    for(int i=0; i<p->knox; ++i)
    {
    xs[p->origin_i+i] += p->XP[i+m];
    xc[p->origin_i+i] += 1.0;
    }
    
    for(int j=0; j<p->knoy; ++j)
    {
    ys[p->origin_j+j] += p->YP[j+m];
    yc[p->origin_j+j] += 1.0;
    }
    
    MPI_Allreduce(MPI_IN_PLACE,xs.data(),gx,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,xc.data(),gx,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,ys.data(),gy,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,yc.data(),gy,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    
    Xg.resize(gx);
    Yg.resize(gy);
    for(int i=0; i<gx; ++i)
    Xg[i] = xs[i]/std::max(xc[i],1.0);
    for(int j=0; j<gy; ++j)
    Yg[j] = ys[j]/std::max(yc[j],1.0);
    
    L.resize(reefcouple::recv_ic.size());
    
    for(size_t n=0; n<L.size(); ++n)
    {
        L[n].ic  = reefcouple::recv_ic[n];
        L[n].app = reefcouple::recv_app[n];
        handshake(p,c,pgc,L[n]);
    }
    
    // bodies driven by a receiver: the first link with I 243
    for(size_t n=0; n<L.size(); ++n)
    if(L[n].nb>0)
    {
        if(body_link<0)
        {
        body_link = int(n);
        pbody = new reefcouple_bodies(this,L[n].nb);
        sixdof_external::attach(pbody);
        
        if(p->mpirank==0)
        cout<<"FNPF coupling: the FNPF bodies 0.."<<L[n].nb-1<<" follow the 6DOF bodies of app "<<L[n].app<<" (I 243)"<<endl;
        
        if(p->X10==0 && p->mpirank==0)
        cout<<"!!! FNPF coupling: I 243 in the receiver, but no FNPF bodies (X 10) !!!"<<endl;
        }
        else
        if(p->mpirank==0)
        cout<<"FNPF coupling: app "<<L[n].app<<" also sends bodies (I 243), only app "<<L[body_link].app<<" drives the FNPF bodies"<<endl;
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
    const int m = increment::marge;
    
    // header to all ranks of the receiver
    const int root = p->mpirank==0 ? MPI_ROOT : MPI_PROC_NULL;
    
    int hi[8] = {reefcouple::version, gx, gy, nz, p->j_dir ? 1 : 0, p->mpi_size, reefcouple::depth, 0};
    MPI_Bcast(hi,8,MPI_INT,root,lk.ic);
    
    std::vector<double> hd(2+gx+gy+nz);
    hd[0] = p->wd;
    hd[1] = p->simtime;
    for(int i=0; i<gx; ++i)
    hd[2+i] = Xg[i];
    for(int j=0; j<gy; ++j)
    hd[2+gx+j] = Yg[j];
    for(int k=0; k<nz; ++k)
    hd[2+gx+gy+k] = p->ZN[k+m];
    
    MPI_Bcast(hd.data(),int(hd.size()),MPI_DOUBLE,root,lk.ic);
    
    // index boxes: mine to the receiver, the receiver windows (and whether it needs Fi) to me
    int nremote;
    MPI_Comm_remote_size(lk.ic,&nremote);
    
    int mine[6] = {i0g,i1g,j0g,j1g,0,0};
    std::vector<int> win(6*nremote);
    MPI_Allgather(mine,6,MPI_INT,win.data(),6,MPI_INT,lk.ic);
    
    lk.fi = (win[4]==1);
    
    lk.d.clear();
    for(int r=0; r<nremote; ++r)
    {
        dest dd;
        dd.rank = r;
        dd.i0 = std::max(i0g,win[6*r]);
        dd.i1 = std::min(i1g,win[6*r+1]);
        dd.j0 = std::max(j0g,win[6*r+2]);
        dd.j1 = std::min(j1g,win[6*r+3]);
        dd.ncol = (dd.i1-dd.i0)*(dd.j1-dd.j0);
        
        if(dd.i1>dd.i0 && dd.j1>dd.j0)
        lk.d.push_back(dd);
    }
    
    const int nf = lk.fi ? 4 : 3;
    const int nd = int(lk.d.size());
    for(int s=0; s<reefcouple::depth; ++s)
    {
        lk.buf[s].resize(nd);
        lk.req[s].assign(nd,MPI_REQUEST_NULL);
        for(int q=0; q<nd; ++q)
        lk.buf[s][q].resize(1+size_t(lk.d[q].ncol)*(1+nf*nz));
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
    
    // feedback configuration from the receiver's rank 0: bodies, overlap columns, levels ahead
    int fc[4];
    MPI_Bcast(fc,4,MPI_INT,0,lk.ic);
    lk.nb   = fc[0];
    lk.nc   = fc[1];
    lk.need = std::max(1,fc[2]);
    lk.fbsize = 2 + 18*lk.nb + 4*lk.nc;
    
    if(lk.nc>0)
    {
        std::vector<double> cl(3*lk.nc);
        MPI_Bcast(cl.data(),3*lk.nc,MPI_DOUBLE,0,lk.ic);
        
        lk.ci.resize(lk.nc);
        lk.cj.resize(lk.nc);
        lk.cr.resize(lk.nc);
        lk.own.clear();
        for(int q=0; q<lk.nc; ++q)
        {
            lk.ci[q] = int(cl[3*q]+0.5);
            lk.cj[q] = int(cl[3*q+1]+0.5);
            lk.cr[q] = cl[3*q+2];
            
            if(lk.ci[q]>=i0g && lk.ci[q]<i1g && lk.cj[q]>=j0g && lk.cj[q]<j1g)
            lk.own.push_back(q);
        }
        lk.te.assign(lk.nc,0.0);
        lk.tf.assign(lk.nc,0.0);
        lk.phin.assign(lk.nc,0.0);
    }
    
    int ntot=0, nlocal=nd;
    MPI_Reduce(&nlocal,&ntot,1,MPI_INT,MPI_SUM,0,pgc->mpi_comm);
    if(p->mpirank==0)
    {
    cout<<"FNPF coupling: app "<<lk.app<<", "<<nremote<<" receiver ranks, "<<ntot<<" FNPF -> receiver messages per step"<<(lk.fi ? ", with Fi (FNPF receiver)" : "")<<endl;
    if(lk.nb>0 || lk.nc>0)
    cout<<"FNPF coupling: app "<<lk.app<<" feeds back "<<lk.nb<<" bodies (I 243) and "<<lk.nc<<" overlap columns (I 244), FNPF waits for it every step"<<endl;
    }
}

void fnpf_couple_send::pack(lexer *p, fdm_fnpf *c, const link &lk, const dest &dd, std::vector<double> &b)
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
        
        if(lk.fi)
        b[1+nc+(3*nc+n)*nz+k] = c->Fi[FIJK];
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
        pack(p,c,lk,lk.d[q],lk.buf[s][q]);
        MPI_Issend(lk.buf[s][q].data(),int(lk.buf[s][q].size()),MPI_DOUBLE,lk.d[q].rank,reefcouple::tag_data,lk.ic,&lk.req[s][q]);
        }
    }
    
    levt.push_back(p->simtime);
    ++level;
}

// rank 0: the messages of the receiver's rank 0 in their order (feedback samples, wait notices, stop)
// until a sample at least as late as tneed has arrived, the receiver waits for a level at a time of at
// least tblock (it cannot send a later sample before the next level), or it has stopped; true: a new
// sample
bool fnpf_couple_send::poll(link &lk, double tneed, double tblock)
{
    bool fresh=false;
    std::vector<double> b;
    
    while(!lk.stopped)
    {
        int flag=0;
        MPI_Status st;
        MPI_Iprobe(0,MPI_ANY_TAG,lk.ic,&flag,&st);
        
        if(!flag)
        {
            const double tl = lk.s2.empty() ? -1.0e20 : lk.s2[0];
            
            if(tl>=tneed-1.0e-10 || lk.twait>=tblock-1.0e-10)
            break;
            
            MPI_Probe(0,MPI_ANY_TAG,lk.ic,&st);
        }
        
        int cnt=0;
        MPI_Get_count(&st,MPI_DOUBLE,&cnt);
        b.resize(std::max(cnt,1));
        MPI_Recv(b.data(),std::max(cnt,1),MPI_DOUBLE,0,st.MPI_TAG,lk.ic,MPI_STATUS_IGNORE);
        
        if(st.MPI_TAG==reefcouple::tag_stop)
        lk.stopped = true;
        
        if(st.MPI_TAG==reefcouple::tag_fbc)
        lk.twait = std::max(lk.twait,b[0]);
        
        if(st.MPI_TAG==reefcouple::tag_fb && cnt==lk.fbsize)
        {
            lk.s1.swap(lk.s2);
            lk.s2 = b;
            fresh = true;
        }
    }
    
    return fresh;
}

void fnpf_couple_send::pre_step(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const int nl = int(L.size());
    const int lev = level-1;          // the last level sent: the current state
    
    // rank 0: feedback samples and stops; flags 0 stop, 1 new sample
    std::vector<int> flag(2*nl,0);
    
    if(p->mpirank==0)
    for(int n=0; n<nl; ++n)
    if(L[n].open)
    {
        link &lk = L[n];
        
        // the receiver needs need levels beyond its time: it has reached at least levt[lev-need]
        double tneed = -1.0e20, tblock = 1.0e30;
        if((lk.nb>0 || lk.nc>0) && lev-lk.need>=0)
        {
        tneed = levt[lev-lk.need];
        tblock = levt[lev-lk.need+1];
        }
        
        flag[2*n+1] = poll(lk,tneed,tblock) ? 1 : 0;
        flag[2*n]   = lk.stopped ? 1 : 0;
    }
    
    MPI_Bcast(flag.data(),2*nl,MPI_INT,0,pgc->mpi_comm);
    
    for(int n=0; n<nl; ++n)
    {
        link &lk = L[n];
        
        if(!lk.open)
        continue;
        
        // the samples to all ranks
        if(flag[2*n+1]==1)
        {
            int sz[2] = {int(lk.s1.size()), int(lk.s2.size())};
            MPI_Bcast(sz,2,MPI_INT,0,pgc->mpi_comm);
            lk.s1.resize(sz[0]);
            lk.s2.resize(sz[1]);
            if(sz[0]>0)
            MPI_Bcast(lk.s1.data(),sz[0],MPI_DOUBLE,0,pgc->mpi_comm);
            if(sz[1]>0)
            MPI_Bcast(lk.s2.data(),sz[1],MPI_DOUBLE,0,pgc->mpi_comm);
        }
        
        if(flag[2*n]==1)
        {
            if(p->mpirank==0)
            cout<<"FNPF coupling: app "<<lk.app<<" has stopped"<<endl;
            
            close(lk);
            lk.targets = false;
            continue;
        }
        
        if(lk.nc>0)
        overlap_targets(p,c,pgc,lk);
    }
}

// targets of the overlap box (I 244) at the end of the step: eta and the surface velocity of the
// receiver extrapolated from its last two samples; Fifsf along the FNPF x lines of the box from
// dPhi/dx = u + w deta/dx, integrated from both ends of a line (the FNPF Fifsf there, extrapolated
// to the end of the step) and blended linearly
void fnpf_couple_send::overlap_targets(lexer *p, fdm_fnpf *c, ghostcell *pgc, link &lk)
{
    const int nc = lk.nc;
    const double tn = p->simtime;
    const double tt = p->simtime + p->dt;
    
    // FNPF Fifsf at the columns now
    std::vector<double> ph(nc,0.0);
    for(int q : lk.own)
    ph[q] = c->Fifsf(lk.ci[q]-p->origin_i,lk.cj[q]-p->origin_j);
    MPI_Allreduce(MPI_IN_PLACE,ph.data(),nc,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    
    std::vector<double> anchor(ph);
    if(lk.tprev>-1.0e19 && tn-lk.tprev>1.0e-12)
    for(int q=0; q<nc; ++q)
    anchor[q] = ph[q] + (ph[q]-lk.phin[q])/(tn-lk.tprev)*(tt-tn);
    
    lk.phin = ph;
    lk.tprev = tn;
    
    if(lk.s2.empty())
    {
    lk.targets = false;
    return;
    }
    
    // receiver data extrapolated to tt
    const int o = 2 + 18*lk.nb;
    const double t2 = lk.s2[0];
    const bool two = (!lk.s1.empty() && t2-lk.s1[0]>1.0e-12);
    const double fac = two ? (tt-t2)/(t2-lk.s1[0]) : 0.0;
    
    std::vector<double> e(nc), u(nc), w(nc);
    for(int q=0; q<nc; ++q)
    {
        e[q] = lk.s2[o+q];
        u[q] = lk.s2[o+nc+q];
        w[q] = lk.s2[o+3*nc+q];
        
        if(two)
        {
        e[q] += fac*(lk.s2[o+q]-lk.s1[o+q]);
        u[q] += fac*(lk.s2[o+nc+q]-lk.s1[o+nc+q]);
        w[q] += fac*(lk.s2[o+3*nc+q]-lk.s1[o+3*nc+q]);
        }
    }
    
    // lines: consecutive columns with the same j and consecutive i
    int q0=0;
    while(q0<nc)
    {
        int q1=q0;
        while(q1+1<nc && lk.cj[q1+1]==lk.cj[q0] && lk.ci[q1+1]==lk.ci[q1]+1)
        ++q1;
        
        const int n = q1-q0+1;
        
        if(n==1)
        lk.tf[q0] = anchor[q0];
        else
        {
            std::vector<double> g(n), fw(n), bw(n);
            for(int a=0; a<n; ++a)
            {
                const int q = q0+a;
                const int qa = (a>0) ? q-1 : q;
                const int qb = (a<n-1) ? q+1 : q;
                const double ex = (e[qb]-e[qa])/(Xg[lk.ci[qb]]-Xg[lk.ci[qa]]);
                g[a] = u[q] + w[q]*ex;
            }
            
            fw[0] = anchor[q0];
            for(int a=1; a<n; ++a)
            fw[a] = fw[a-1] + 0.5*(g[a-1]+g[a])*(Xg[lk.ci[q0+a]]-Xg[lk.ci[q0+a-1]]);
            
            bw[n-1] = anchor[q1];
            for(int a=n-2; a>=0; --a)
            bw[a] = bw[a+1] - 0.5*(g[a+1]+g[a])*(Xg[lk.ci[q0+a+1]]-Xg[lk.ci[q0+a]]);
            
            const double xa = Xg[lk.ci[q0]], xb = Xg[lk.ci[q1]];
            for(int a=0; a<n; ++a)
            {
            const double s = (Xg[lk.ci[q0+a]]-xa)/(xb-xa);
            lk.tf[q0+a] = (1.0-s)*fw[a] + s*bw[a];
            }
        }
        
        for(int q=q0; q<=q1; ++q)
        lk.te[q] = e[q];
        
        q0 = q1+1;
    }
    
    lk.targets = true;
}

void fnpf_couple_send::relax(lexer *p, slice &f, bool eta)
{
    for(link &lk : L)
    if(lk.open && lk.targets)
    for(int q : lk.own)
    {
        const int i = lk.ci[q]-p->origin_i;
        const int j = lk.cj[q]-p->origin_j;
        const double r = lk.cr[q];
        const double tg = eta ? lk.te[q] : lk.tf[q];
        
        f(i,j) = r*tg + (1.0-r)*f(i,j);
    }
}

void fnpf_couple_send::relax_eta(lexer *p, slice &f)
{
    if(reefcouple::sender!=nullptr)
    reefcouple::sender->relax(p,f,true);
}

void fnpf_couple_send::relax_fifsf(lexer *p, slice &f)
{
    if(reefcouple::sender!=nullptr)
    reefcouple::sender->relax(p,f,false);
}

bool fnpf_couple_send::body_state(int n, double t, double *s) const
{
    if(body_link<0)
    return false;
    
    const link &lk = L[body_link];
    
    if(lk.s2.empty() || n>=lk.nb)
    return false;
    
    const double *b2 = &lk.s2[2+18*n];
    const double dl = t-lk.s2[0];
    
    // acceleration and angular acceleration from the last two samples
    double a[3]={0.0,0.0,0.0}, al[3]={0.0,0.0,0.0};
    if(!lk.s1.empty() && lk.s2[0]-lk.s1[0]>1.0e-12)
    {
        const double *b1 = &lk.s1[2+18*n];
        const double h = lk.s2[0]-lk.s1[0];
        for(int q=0; q<3; ++q)
        {
        a[q]  = (b2[3+q]-b1[3+q])/h;
        al[q] = (b2[6+q]-b1[6+q])/h;
        }
    }
    
    for(int q=0; q<3; ++q)
    {
    s[q]   = b2[q] + b2[3+q]*dl + 0.5*a[q]*dl*dl;
    s[3+q] = b2[3+q] + a[q]*dl;
    s[6+q] = b2[6+q] + al[q]*dl;
    }
    
    // rotation: R(t) = exp([th]) R2, th = (w2 + al dl/2) dl (inertial frame, Rodrigues)
    Eigen::Vector3d th;
    for(int q=0; q<3; ++q)
    th(q) = (b2[6+q] + 0.5*al[q]*dl)*dl;
    
    Eigen::Matrix3d R2;
    for(int r=0; r<3; ++r)
    for(int q=0; q<3; ++q)
    R2(r,q) = b2[9+3*r+q];
    
    Eigen::Matrix3d Rd = Eigen::Matrix3d::Identity();
    const double ang = th.norm();
    if(ang>1.0e-14)
    {
        const Eigen::Vector3d k = th/ang;
        Eigen::Matrix3d K;
        K <<    0.0, -k(2),  k(1),
               k(2),   0.0, -k(0),
              -k(1),  k(0),   0.0;
        Rd = Eigen::Matrix3d::Identity() + sin(ang)*K + (1.0-cos(ang))*K*K;
    }
    
    const Eigen::Matrix3d R = Rd*R2;
    for(int r=0; r<3; ++r)
    for(int q=0; q<3; ++q)
    s[9+3*r+q] = R(r,q);
    
    return true;
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
    {
        if(lk.open)
        close(lk);
        
        // the receiver sends its stop when it ends: its remaining samples, then the stop
        if(p->mpirank==0)
        while(!lk.stopped)
        poll(lk,1.0e30,1.0e30);
    }
    
    if(pbody!=nullptr && sixdof_external::instance()==pbody)
    sixdof_external::attach(nullptr);
    
    finished = true;
}
