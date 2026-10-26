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

#include"wave_lib_fnpf_couple.h"
#include"reefcouple.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iostream>
#include<algorithm>
#include<cmath>

wave_lib_fnpf_couple::wave_lib_fnpf_couple(lexer *pp, ghostcell *pgc) : wave_lib_parameters(pp,pgc), p(pp)
{
    ic = reefcouple::src_ic;
    
    if(ic==MPI_COMM_NULL)
    {
        if(p->mpirank==0)
        cout<<"!!! B 92 62: no FNPF source in this MPI job !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    if(reefcouple::recv!=nullptr)
    {
        if(p->mpirank==0)
        cout<<"!!! B 92 62: one FNPF coupling per model !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    // the decomposed (space/time) evaluation is for wave theories only
    p->B89 = 0;
    
    order = (p->I242==1) ? 1 : 3;
    
    // frame: wave frame of B 105 (as for the HDC, origin 0,0 unless B 105 is given), then the offsets
    const double g = p->B105_1*(PI/180.0);
    cg = cos(g);
    sg = sin(g);
    dx0 = p->I231;
    dy0 = p->I232;
    dz0 = p->I233;
    
    handshake(p,pgc);
    
    reefcouple::recv = this;
    
    // the state at the start
    tlast = p->simtime + p->I241;
    weights(tlast);
}

wave_lib_fnpf_couple::~wave_lib_fnpf_couple()
{
}

void wave_lib_fnpf_couple::handshake(lexer *p, ghostcell *pgc)
{
    // FNPF grid
    int hi[8];
    MPI_Bcast(hi,8,MPI_INT,0,ic);
    
    if(hi[0]!=reefcouple::version)
    {
        if(p->mpirank==0)
        cout<<"!!! B 92 62: coupling protocol "<<hi[0]<<" of the FNPF model, "<<reefcouple::version<<" here !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    gx    = hi[1];
    gy    = hi[2];
    nz    = hi[3];
    jdirF = hi[4];
    nfnpf = hi[5];
    
    std::vector<double> hd(2+gx+gy+nz);
    MPI_Bcast(hd.data(),int(hd.size()),MPI_DOUBLE,0,ic);
    
    wdF = hd[0];
    t0F = hd[1];
    X.assign(hd.begin()+2,hd.begin()+2+gx);
    Y.assign(hd.begin()+2+gx,hd.begin()+2+gx+gy);
    ZN.assign(hd.begin()+2+gx+gy,hd.end());
    
    // window: the subdomain with a margin, in the FNPF frame
    const double e = 4.0*p->DXM;
    const double xa = p->originx-e, xb = p->endx+e;
    const double ya = p->originy-e, yb = p->endy+e;
    const double cx[4] = {xa,xb,xa,xb};
    const double cy[4] = {ya,ya,yb,yb};
    
    double xmin=1.0e20, xmax=-1.0e20, ymin=1.0e20, ymax=-1.0e20;
    for(int q=0; q<4; ++q)
    {
        const double xr = cx[q]-p->B105_2;
        const double yr = cy[q]-p->B105_3;
        const double xf =  xr*cg + yr*sg + dx0;
        const double yf = -xr*sg + yr*cg + dy0;
        xmin = std::min(xmin,xf);
        xmax = std::max(xmax,xf);
        ymin = std::min(ymin,yf);
        ymax = std::max(ymax,yf);
    }
    
    auto range = [](const std::vector<double> &A, double lo, double hi, int &a, int &b)
    {
        const int n = int(A.size());
        a = int(std::upper_bound(A.begin(),A.end(),lo)-A.begin())-2;
        b = int(std::lower_bound(A.begin(),A.end(),hi)-A.begin())+2;
        a = std::max(0,std::min(a,n-1));
        b = std::max(a+1,std::min(b,n));
    };
    
    range(X,xmin,xmax,wi0,wi1);
    
    if(jdirF==1)
    range(Y,ymin,ymax,wj0,wj1);
    else
    {
    wj0 = 0;
    wj1 = gy;
    }
    
    nwx = wi1-wi0;
    nwy = wj1-wj0;
    nwc = nwx*nwy;
    
    // FNPF boxes <-> window
    int mine[4] = {wi0,wi1,wj0,wj1};
    std::vector<int> box(4*nfnpf);
    MPI_Allgather(mine,4,MPI_INT,box.data(),4,MPI_INT,ic);
    
    S.clear();
    for(int r=0; r<nfnpf; ++r)
    {
        source s;
        s.rank = r;
        s.i0 = std::max(box[4*r],wi0);
        s.i1 = std::min(box[4*r+1],wi1);
        s.j0 = std::max(box[4*r+2],wj0);
        s.j1 = std::min(box[4*r+3],wj1);
        s.ncol = (s.i1-s.i0)*(s.j1-s.j0);
        
        if(s.i1>s.i0 && s.j1>s.j0)
        {
        s.buf.resize(1+size_t(s.ncol)*(1+3*nz));
        S.push_back(s);
        }
    }
    
    // bed
    B.assign(nwc,0.0);
    {
        std::vector<std::vector<double>> bb(S.size());
        std::vector<MPI_Request> rq(S.size());
        for(size_t q=0; q<S.size(); ++q)
        {
        bb[q].resize(S[q].ncol);
        MPI_Irecv(bb[q].data(),S[q].ncol,MPI_DOUBLE,S[q].rank,reefcouple::tag_bed,ic,&rq[q]);
        }
        MPI_Waitall(int(rq.size()),rq.data(),MPI_STATUSES_IGNORE);
        
        for(size_t q=0; q<S.size(); ++q)
        {
            int n=0;
            for(int ig=S[q].i0; ig<S[q].i1; ++ig)
            for(int jg=S[q].j0; jg<S[q].j1; ++jg)
            B[(ig-wi0)*nwy + (jg-wj0)] = bb[q][n++];
        }
    }
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: FNPF coupling (B 92 62)"<<endl;
    cout<<" FNPF grid: "<<gx<<" x "<<gy<<" x "<<nz<<" nodes, "<<(jdirF==1 ? "3D" : "2D")<<", wd: "<<wdF<<", t0: "<<t0F<<endl;
    cout<<" receiver: "<<(p->j_dir==1 ? "3D" : "2D")<<", B 105 angle: "<<p->B105_1<<" origin: "<<p->B105_2<<" "<<p->B105_3
        <<", offsets I 231-233: "<<dx0<<" "<<dy0<<" "<<dz0<<", I 241: "<<p->I241<<endl;
    cout<<" time interpolation: "<<(order==3 ? "cubic" : "linear")<<endl;
    }
    
    if(p->j_dir==1 && jdirF==0 && p->mpirank==0)
    cout<<" 2D FNPF -> 3D receiver: long-crested waves in the B 105 direction"<<endl;
}

bool wave_lib_fnpf_couple::recv_level()
{
    if(ended)
    return false;
    
    std::vector<MPI_Request> rq(S.size());
    std::vector<MPI_Status> st(S.size());
    for(size_t q=0; q<S.size(); ++q)
    MPI_Irecv(S[q].buf.data(),int(S[q].buf.size()),MPI_DOUBLE,S[q].rank,MPI_ANY_TAG,ic,&rq[q]);
    
    MPI_Waitall(int(rq.size()),rq.data(),st.data());
    
    bool any_end=false;
    for(size_t q=0; q<S.size(); ++q)
    if(st[q].MPI_TAG==reefcouple::tag_end)
    {
    S[q].ended = true;
    any_end = true;
    }
    
    // the FNPF model has ended: the last level is held
    if(any_end)
    {
        ended = true;
        if(p->mpirank==0)
        cout<<"Wave_Lib: the FNPF model has ended, its last state is held"<<endl;
        return false;
    }
    
    level L;
    L.t = S[0].buf[0];
    L.E.assign(nwc,0.0);
    L.U.assign(size_t(nwc)*nz,0.0);
    L.V.assign(size_t(nwc)*nz,0.0);
    L.W.assign(size_t(nwc)*nz,0.0);
    
    for(source &s : S)
    {
        const size_t nc = s.ncol;
        const double *b = s.buf.data();
        size_t n=0;
        
        for(int ig=s.i0; ig<s.i1; ++ig)
        for(int jg=s.j0; jg<s.j1; ++jg)
        {
            const size_t col = size_t(ig-wi0)*nwy + (jg-wj0);
            
            L.E[col] = b[1+n];
            
            for(int k=0; k<nz; ++k)
            {
            L.U[col*nz+k] = b[1+nc+n*nz+k];
            L.V[col*nz+k] = b[1+nc+(nc+n)*nz+k];
            L.W[col*nz+k] = b[1+nc+(2*nc+n)*nz+k];
            }
            ++n;
        }
    }
    
    lev.push_back(std::move(L));
    
    return true;
}

void wave_lib_fnpf_couple::weights(double t)
{
    const int need = (order==3) ? 2 : 1;
    
    auto ahead = [&]()
    {
        int n=0;
        for(const level &l : lev)
        if(l.t>t)
        ++n;
        return n;
    };
    
    // the levels up to t (and one more for cubic)
    while(!ended && ahead()<need)
    if(!recv_level())
    break;
    
    // drop the levels no longer needed (one spare for a slightly earlier request)
    int behind = int(lev.size())-ahead();
    while(behind>need+1 && lev.size()>1)
    {
    lev.pop_front();
    --behind;
    }
    
    wl.clear();
    ww_.clear();
    
    const int n = int(lev.size());
    
    if(n==0)
    return;
    
    int q=-1;
    for(int m=0; m<n; ++m)
    if(lev[m].t<=t)
    q=m;
    
    // before the first level or after the last (FNPF ended): hold
    if(q<0 || q==n-1)
    {
    wl.push_back(q<0 ? 0 : n-1);
    ww_.push_back(1.0);
    return;
    }
    
    // cubic Lagrange on q-1 .. q+2
    if(order==3 && q>=1 && q+2<=n-1)
    {
        double tt[4];
        for(int a=0; a<4; ++a)
        tt[a] = lev[q-1+a].t;
        
        for(int a=0; a<4; ++a)
        {
            double w=1.0;
            for(int b=0; b<4; ++b)
            if(b!=a)
            w *= (t-tt[b])/(tt[a]-tt[b]);
            
            wl.push_back(q-1+a);
            ww_.push_back(w);
        }
        return;
    }
    
    // linear
    const double dt = lev[q+1].t-lev[q].t;
    const double a = dt>1.0e-20 ? (t-lev[q].t)/dt : 1.0;
    wl.push_back(q);
    ww_.push_back(1.0-a);
    wl.push_back(q+1);
    ww_.push_back(a);
}

wave_lib_fnpf_couple::hstencil wave_lib_fnpf_couple::locate(double xf, double yf) const
{
    hstencil h;
    int ia=wi0, ja=wj0;
    double ax=0.0, ay=0.0;
    
    if(nwx>1)
    {
        ia = int(std::upper_bound(X.begin()+wi0,X.begin()+wi1,xf)-X.begin())-1;
        
        if(ia<wi0)
        {
        ia = wi0;
        ax = 0.0;
        }
        else
        if(ia>=wi1-1)
        {
        ia = wi1-2;
        ax = 1.0;
        }
        else
        ax = (xf-X[ia])/(X[ia+1]-X[ia]);
    }
    
    if(jdirF==1 && nwy>1)
    {
        ja = int(std::upper_bound(Y.begin()+wj0,Y.begin()+wj1,yf)-Y.begin())-1;
        
        if(ja<wj0)
        {
        ja = wj0;
        ay = 0.0;
        }
        else
        if(ja>=wj1-1)
        {
        ja = wj1-2;
        ay = 1.0;
        }
        else
        ay = (yf-Y[ja])/(Y[ja+1]-Y[ja]);
    }
    
    const int ib = (nwx>1) ? ia+1 : ia;
    const int jb = (jdirF==1 && nwy>1) ? ja+1 : ja;
    
    h.col[0] = (ia-wi0)*nwy + (ja-wj0);
    h.col[1] = (ib-wi0)*nwy + (ja-wj0);
    h.col[2] = (ia-wi0)*nwy + (jb-wj0);
    h.col[3] = (ib-wi0)*nwy + (jb-wj0);
    
    h.w[0] = (1.0-ax)*(1.0-ay);
    h.w[1] = ax*(1.0-ay);
    h.w[2] = (1.0-ax)*ay;
    h.w[3] = ax*ay;
    
    return h;
}

double wave_lib_fnpf_couple::eta_h(const hstencil &h) const
{
    double val=0.0;
    
    for(int c=0; c<4; ++c)
    if(h.w[c]>0.0)
    {
        double e=0.0;
        for(size_t m=0; m<wl.size(); ++m)
        e += ww_[m]*lev[wl[m]].E[h.col[c]];
        
        val += h.w[c]*e;
    }
    
    return val;
}

void wave_lib_fnpf_couple::vel_h(const hstencil &h, double zf, double &u, double &v, double &w) const
{
    u=v=w=0.0;
    
    for(int c=0; c<4; ++c)
    if(h.w[c]>0.0)
    {
        const int col = h.col[c];
        
        double e=0.0;
        for(size_t m=0; m<wl.size(); ++m)
        e += ww_[m]*lev[wl[m]].E[col];
        
        const double depth = e + wdF - B[col];
        
        if(depth<1.0e-6)
        continue;
        
        // sigma of zf in this column, the nodes k, k+1 around it
        double s = (zf-B[col])/depth;
        s = std::max(ZN.front(),std::min(s,ZN.back()));
        
        int k=0;
        double b=0.0;
        if(nz>1)
        {
            k = int(std::upper_bound(ZN.begin(),ZN.end(),s)-ZN.begin())-1;
            k = std::max(0,std::min(k,nz-2));
            b = (s-ZN[k])/(ZN[k+1]-ZN[k]);
        }
        const int kp = (nz>1) ? k+1 : k;
        
        double uc=0.0, vc=0.0, wc=0.0;
        for(size_t m=0; m<wl.size(); ++m)
        {
            const level &L = lev[wl[m]];
            const size_t n0 = size_t(col)*nz + k;
            const size_t n1 = size_t(col)*nz + kp;
            uc += ww_[m]*(L.U[n0]*(1.0-b) + L.U[n1]*b);
            vc += ww_[m]*(L.V[n0]*(1.0-b) + L.V[n1]*b);
            wc += ww_[m]*(L.W[n0]*(1.0-b) + L.W[n1]*b);
        }
        
        u += h.w[c]*uc;
        v += h.w[c]*vc;
        w += h.w[c]*wc;
    }
    
    // FNPF 2D: no transverse velocity
    if(jdirF==0)
    v=0.0;
}

void wave_lib_fnpf_couple::rotate(double uf, double vf, double &u, double &v) const
{
    u = cg*uf - sg*vf;
    v = sg*uf + cg*vf;
}

void wave_lib_fnpf_couple::frame(double x, double y, double &xf, double &yf) const
{
    // x, y: the wave frame of B 105 (iowave xgen, ygen)
    xf = x + dx0;
    yf = y + dy0;
}

double wave_lib_fnpf_couple::wave_eta(lexer *p, double x, double y)
{
    sync(p);
    double xf,yf;
    frame(x,y,xf,yf);
    return eta_h(locate(xf,yf));
}

double wave_lib_fnpf_couple::wave_u(lexer *p, double x, double y, double z)
{
    sync(p);
    double xf,yf,uf,vf,wf,u,v;
    frame(x,y,xf,yf);
    vel_h(locate(xf,yf),z+wdF+dz0,uf,vf,wf);
    rotate(uf,vf,u,v);
    return u;
}

double wave_lib_fnpf_couple::wave_v(lexer *p, double x, double y, double z)
{
    sync(p);
    double xf,yf,uf,vf,wf,u,v;
    frame(x,y,xf,yf);
    vel_h(locate(xf,yf),z+wdF+dz0,uf,vf,wf);
    rotate(uf,vf,u,v);
    return v;
}

double wave_lib_fnpf_couple::wave_w(lexer *p, double x, double y, double z)
{
    sync(p);
    double xf,yf,uf,vf,wf;
    frame(x,y,xf,yf);
    vel_h(locate(xf,yf),z+wdF+dz0,uf,vf,wf);
    return wf;
}

double wave_lib_fnpf_couple::wave_fi(lexer *p, double x, double y, double z)
{
    return 0.0;
}

void wave_lib_fnpf_couple::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_fnpf_couple::wave_prestep(lexer *p, ghostcell *pgc)
{
    tlast = p->wavetime + p->I241;
    weights(tlast);
}

// a query at another wave time than the last (sub-steps, e.g. CFD AMR patches): new time weights
void wave_lib_fnpf_couple::sync(lexer *p)
{
    const double t = p->wavetime + p->I241;
    
    if(t!=tlast)
    {
    tlast = t;
    weights(t);
    }
}

void wave_lib_fnpf_couple::wave_cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    cache_x = x;
    cache_y = y;
    
    cache.resize(x.size());
    for(size_t q=0; q<x.size(); ++q)
    {
        double xf,yf;
        frame(x[q],y[q],xf,yf);
        cache[q] = locate(xf,yf);
    }
}

double wave_lib_fnpf_couple::wave_eta_c(lexer *p, int q)
{
    sync(p);
    return eta_h(cache[q]);
}

double wave_lib_fnpf_couple::wave_u_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return u;
}

double wave_lib_fnpf_couple::wave_v_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return v;
}

double wave_lib_fnpf_couple::wave_w_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return w;
}

void wave_lib_fnpf_couple::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    sync(p);
    double uf,vf;
    vel_h(cache[q],z+wdF+dz0,uf,vf,w);
    rotate(uf,vf,u,v);
}

void wave_lib_fnpf_couple::finish()
{
    if(finished)
    return;
    
    finished = true;
    
    // stop to the source (rank 0 for the model), then everything the source still sends
    if(p->mpirank==0)
    {
    int v=1;
    MPI_Send(&v,1,MPI_INT,0,reefcouple::tag_stop,ic);
    }
    
    int open=0;
    for(const source &s : S)
    if(!s.ended)
    ++open;
    
    std::vector<double> tmp;
    while(open>0)
    {
        MPI_Status st;
        MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,ic,&st);
        int cnt=0;
        MPI_Get_count(&st,MPI_DOUBLE,&cnt);
        tmp.resize(std::max(cnt,1));
        MPI_Recv(tmp.data(),cnt,MPI_DOUBLE,st.MPI_SOURCE,st.MPI_TAG,ic,MPI_STATUS_IGNORE);
        
        if(st.MPI_TAG==reefcouple::tag_end)
        for(source &s : S)
        if(s.rank==st.MPI_SOURCE && !s.ended)
        {
        s.ended = true;
        --open;
        }
    }
}
