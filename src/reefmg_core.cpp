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

#include "reefmg_core.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>

reefmg_core::reefmg_core()
{
    comm=MPI_COMM_NULL;
    myrank=nprocs=0;
    nbx0=nbx1=nby0=nby1=MPI_PROC_NULL;
    coarse_sweeps=16;
    nfallback=0;
    sweepstyle=0;
    errmsg[0]='\0';
}

reefmg_core::~reefmg_core()
{
    //  the solver object may outlive MPI if it is torn down late
    int fin=0;
    MPI_Finalized(&fin);

    if(comm!=MPI_COMM_NULL && fin==0)
    MPI_Comm_free(&comm);
}

bool reefmg_core::setup(MPI_Comm world,int npx,int npy,int cx,int cy,
                             int nx,int ny,int nz,int gnx,int gny,int maxlevel,
                             const double *dxn,const double *dyn)
{
    //  Re-order the ranks so that neighbours are a simple offset apart.
    const int key=cy*npx+cx;
    MPI_Comm_split(world,0,key,&comm);
    MPI_Comm_rank(comm,&myrank);
    MPI_Comm_size(comm,&nprocs);

    if(myrank!=key)
    {
        snprintf(errmsg,sizeof(errmsg),"rank re-ordering failed");
        return false;
    }

    nbx0 = (cx>0    )? myrank-1   : MPI_PROC_NULL;
    nbx1 = (cx<npx-1)? myrank+1   : MPI_PROC_NULL;
    nby0 = (cy>0    )? myrank-npx : MPI_PROC_NULL;
    nby1 = (cy<npy-1)? myrank+npx : MPI_PROC_NULL;

    //  How far can every rank coarsen?  All ranks must agree, so take the
    //  global minimum.  Coarsening stops when a local extent turns odd or
    //  drops below 8 cells.
    //  x and y are coarsened independently, so a 2D vertical run (ny==1)
    //  simply keeps y at one cell on every level.
    int mlx=0,mly=0;
    {
        int ax=nx; while(ax>=8 && (maxlevel<=0 || mlx<maxlevel)){ax=(ax+1)/2; ++mlx;}
        int ay=ny; while(ay>=8 && (maxlevel<=0 || mly<maxlevel)){ay=(ay+1)/2; ++mly;}
    }
    int levx=0,levy=0;
    MPI_Allreduce(&mlx,&levx,1,MPI_INT,MPI_MIN,comm);
    MPI_Allreduce(&mly,&levy,1,MPI_INT,MPI_MIN,comm);
    const int nlev=std::max(levx,levy);

    lev.resize(nlev+1);
    int ax=nx, ay=ny, agx=gnx, agy=gny;

    for(int l=0;l<=nlev;++l)
    {
        sc_level &L=lev[l];
        L.lid=l;
        L.nx=ax; L.ny=ay; L.nz=nz; L.gnx=agx; L.gny=agy;
        L.rx=(l>0 && l<=levx)?2:1;
        L.ry=(l>0 && l<=levy)?2:1;
        const long N=L.size();
        L.p.assign(N,0.0); L.n.assign(N,0.0); L.s.assign(N,0.0);
        L.w.assign(N,0.0); L.e.assign(N,0.0); L.t.assign(N,0.0); L.b.assign(N,0.0);
        L.u.assign(N,0.0); L.f.assign(N,0.0); L.r.assign(N,0.0);
        L.act.assign(N,0);
        if(l<levx){ax=(ax+1)/2; agx=(agx+1)/2;}
        if(l<levy){ay=(ay+1)/2; agy=(agy+1)/2;}
    }

    const long fn=lev[0].size();
    kr.assign(fn,0.0); krhat.assign(fn,0.0); kp.assign(fn,0.0); kv.assign(fn,0.0);
    ks.assign(fn,0.0); kt.assign(fn,0.0);    ky.assign(fn,0.0); kz.assign(fn,0.0);

    const long mx=std::max((long)lev[0].ny,(long)lev[0].nx)*nz;
    sbuf.assign(2*mx,0.0); rbuf.assign(2*mx,0.0);

    build_widths(dxn,dyn);

    return true;
}

//  Cell widths per level.  A coarse cell is as wide as its children together,
//  which is what makes a ragged agglomerate - one child instead of two - come
//  out with the right volume rather than a quarter of it.
void reefmg_core::build_widths(const double *dxn,const double *dyn)
{
    sc_level &F=lev[0];
    F.hx.assign(F.nx+2,1.0);
    F.hy.assign(F.ny+2,1.0);

    if(dxn) for(int i=0;i<F.nx+2;++i) F.hx[i]=dxn[i];
    if(dyn) for(int j=0;j<F.ny+2;++j) F.hy[j]=dyn[j];

    for(int l=0;l+1<(int)lev.size();++l)
    {
        sc_level &C=lev[l+1];
        sc_level &P=lev[l];

        C.hx.assign(C.nx+2,0.0);
        C.hy.assign(C.ny+2,0.0);

        for(int I=0;I<C.nx;++I)
        {
            double s=0.0;
            for(int a=0;a<C.rx;++a)
            {
                const int i=C.rx*I+a;
                if(i<P.nx) s+=P.hx[i+1];
            }
            C.hx[I+1]=(s>0.0? s : 1.0);
        }
        for(int J=0;J<C.ny;++J)
        {
            double s=0.0;
            for(int b=0;b<C.ry;++b)
            {
                const int j=C.ry*J+b;
                if(j<P.ny) s+=P.hy[j+1];
            }
            C.hy[J+1]=(s>0.0? s : 1.0);
        }

        exchange_widths(C);
    }

    exchange_widths(F);
}

//  one cell of width halo, so the coarse centre distance at a process
//  boundary is the real one
void reefmg_core::exchange_widths(sc_level &L)
{
    double sx[2]={L.hx[1],L.hx[L.nx]}, rx[2]={L.hx[1],L.hx[L.nx]};
    MPI_Sendrecv(&sx[0],1,MPI_DOUBLE,nbx0,701,&rx[1],1,MPI_DOUBLE,nbx1,701,comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&sx[1],1,MPI_DOUBLE,nbx1,702,&rx[0],1,MPI_DOUBLE,nbx0,702,comm,MPI_STATUS_IGNORE);
    L.hx[0]=rx[0]; L.hx[L.nx+1]=rx[1];

    double sy[2]={L.hy[1],L.hy[L.ny]}, ry[2]={L.hy[1],L.hy[L.ny]};
    MPI_Sendrecv(&sy[0],1,MPI_DOUBLE,nby0,703,&ry[1],1,MPI_DOUBLE,nby1,703,comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&sy[1],1,MPI_DOUBLE,nby1,704,&ry[0],1,MPI_DOUBLE,nby0,704,comm,MPI_STATUS_IGNORE);
    L.hy[0]=ry[0]; L.hy[L.ny+1]=ry[1];
}

//  ------------------------------------------------------------------ halo
//  A plane of constant i is contiguous in memory, a plane of constant j is
//  strided, so the second one is packed.
void reefmg_core::halo(sc_level &L)
{
    if(nbx0==MPI_PROC_NULL && nbx1==MPI_PROC_NULL &&
       nby0==MPI_PROC_NULL && nby1==MPI_PROC_NULL)
    return;

    //  tags kept well apart per level so that nothing can be confused even if
    //  a future change overlaps two exchanges
    const int tx0=1000+4*L.lid, tx1=tx0+1, ty0=5000+4*L.lid, ty1=ty0+1;

    MPI_Request req[8]; int nreq=0;
    const int nz=L.nz;

    //  x direction: contiguous blocks of ny*nz
    const int cntx=L.ny*nz;
    double *sx0=&L.u[L.idx(0,0,0)];
    double *sx1=&L.u[L.idx(L.nx-1,0,0)];
    double *rx0=&L.u[L.idx(-1,0,0)];
    double *rx1=&L.u[L.idx(L.nx,0,0)];

    MPI_Irecv(rx0,cntx,MPI_DOUBLE,nbx0,tx0,comm,&req[nreq++]);
    MPI_Irecv(rx1,cntx,MPI_DOUBLE,nbx1,tx1,comm,&req[nreq++]);
    MPI_Isend(sx0,cntx,MPI_DOUBLE,nbx0,tx1,comm,&req[nreq++]);
    MPI_Isend(sx1,cntx,MPI_DOUBLE,nbx1,tx0,comm,&req[nreq++]);
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    //  y direction: pack
    nreq=0;
    const int cnty=L.nx*nz;
    double *p0=&sbuf[0], *p1=&sbuf[cnty], *q0=&rbuf[0], *q1=&rbuf[cnty];

    for(int i=0;i<L.nx;++i) for(int k=0;k<nz;++k)
    {
        p0[i*nz+k]=L.u[L.idx(i,0,k)];
        p1[i*nz+k]=L.u[L.idx(i,L.ny-1,k)];
    }

    MPI_Irecv(q0,cnty,MPI_DOUBLE,nby0,ty0,comm,&req[nreq++]);
    MPI_Irecv(q1,cnty,MPI_DOUBLE,nby1,ty1,comm,&req[nreq++]);
    MPI_Isend(p0,cnty,MPI_DOUBLE,nby0,ty1,comm,&req[nreq++]);
    MPI_Isend(p1,cnty,MPI_DOUBLE,nby1,ty0,comm,&req[nreq++]);
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    if(nby0!=MPI_PROC_NULL)
    for(int i=0;i<L.nx;++i) for(int k=0;k<nz;++k) L.u[L.idx(i,-1,k)]=q0[i*nz+k];
    if(nby1!=MPI_PROC_NULL)
    for(int i=0;i<L.nx;++i) for(int k=0;k<nz;++k) L.u[L.idx(i,L.ny,k)]=q1[i*nz+k];
}

//  --------------------------------------------------------- coarse operators
//  Horizontal: a coarse face is two fine faces, so the coarse coefficient is
//  an eighth of their sum (half the conductance over four times the volume).
//  Vertical: a coarse face is four fine faces over four times the volume, so
//  the coefficient is their mean.  The row sum is carried across unchanged,
//  which is what keeps the free-surface Dirichlet term alive on every level.
void reefmg_core::coarsen()
{
    for(int l=0;l+1<(int)lev.size();++l)
    {
        sc_level &F=lev[l];
        sc_level &C=lev[l+1];

        std::fill(C.p.begin(),C.p.end(),0.0);
        std::fill(C.n.begin(),C.n.end(),0.0); std::fill(C.s.begin(),C.s.end(),0.0);
        std::fill(C.w.begin(),C.w.end(),0.0); std::fill(C.e.begin(),C.e.end(),0.0);
        std::fill(C.t.begin(),C.t.end(),0.0); std::fill(C.b.begin(),C.b.end(),0.0);
        std::fill(C.act.begin(),C.act.end(),0);

        for(int I=0;I<C.nx;++I)
        for(int J=0;J<C.ny;++J)
        for(int k=0;k<C.nz;++k)
        {
            const long qc=C.idx(I,J,k);

            const int rx=C.rx, ry=C.ry;

            //  Recover the face coefficient a from the fine row
            //     M = -a/(d*h)      d = centre distance, h = cell width
            //  average it over the fine faces that make up the coarse face,
            //  then rebuild with the coarse width and centre distance.  A
            //  ragged agglomerate simply has a smaller H, and a stretched
            //  grid a varying one; neither needs a special case.
            double volsum=0.0, tsum=0.0, bsum=0.0, xsum=0.0;
            double an=0.0,as=0.0,aw=0.0,ae=0.0;
            double wn=0.0,ws=0.0,ww_=0.0,we=0.0;
            int na=0, ihi=rx*I, jhi=ry*J;

            for(int a=0;a<rx;++a) if(rx*I+a<F.nx) ihi=rx*I+a;
            for(int b=0;b<ry;++b) if(ry*J+b<F.ny) jhi=ry*J+b;

            for(int i=rx*I;i<=ihi;++i)
            for(int j=ry*J;j<=jhi;++j)
            {
                const long qf=F.idx(i,j,k);
                if(F.act[qf]==0) continue;
                ++na;

                const double hxf=F.hx[i+1], hyf=F.hy[j+1];
                const double vol=hxf*hyf;

                volsum+=vol;
                tsum  +=F.t[qf]*vol;
                bsum  +=F.b[qf]*vol;
                xsum  +=(F.p[qf]+F.n[qf]+F.s[qf]+F.w[qf]+F.e[qf]+F.t[qf]+F.b[qf])*vol;

                if(i==ihi)
                {
                    const double d=0.5*(hxf+F.hx[i+2]);
                    an+=(-F.n[qf]*d*hxf)*hyf; wn+=hyf;
                }
                if(i==rx*I)
                {
                    const double d=0.5*(hxf+F.hx[i]);
                    as+=(-F.s[qf]*d*hxf)*hyf; ws+=hyf;
                }
                if(j==jhi)
                {
                    const double d=0.5*(hyf+F.hy[j+2]);
                    aw+=(-F.w[qf]*d*hyf)*hxf; ww_+=hxf;
                }
                if(j==ry*J)
                {
                    const double d=0.5*(hyf+F.hy[j]);
                    ae+=(-F.e[qf]*d*hyf)*hxf; we+=hxf;
                }
            }

            if(na==0)
            {
                C.p[qc]=1.0;
                continue;
            }

            C.act[qc]=1;

            const double HX=C.hx[I+1], HY=C.hy[J+1];
            const double DXn=0.5*(HX+C.hx[I+2]), DXs=0.5*(HX+C.hx[I]);
            const double DYw=0.5*(HY+C.hy[J+2]), DYe=0.5*(HY+C.hy[J]);

            C.t[qc]=tsum/volsum;
            C.b[qc]=bsum/volsum;
            C.n[qc]=(wn >0.0)? -(an /wn )/(DXn*HX) : 0.0;
            C.s[qc]=(ws >0.0)? -(as /ws )/(DXs*HX) : 0.0;
            C.w[qc]=(ww_>0.0)? -(aw /ww_)/(DYw*HY) : 0.0;
            C.e[qc]=(we >0.0)? -(ae /we )/(DYe*HY) : 0.0;

            if(I==0        && nbx0==MPI_PROC_NULL) C.s[qc]=0.0;
            if(I==C.nx-1   && nbx1==MPI_PROC_NULL) C.n[qc]=0.0;
            if(J==0        && nby0==MPI_PROC_NULL) C.e[qc]=0.0;
            if(J==C.ny-1   && nby1==MPI_PROC_NULL) C.w[qc]=0.0;
            if(k==0)        C.b[qc]=0.0;
            if(k==C.nz-1)   C.t[qc]=0.0;

            const double excess=xsum/volsum;
            C.p[qc]=excess-(C.n[qc]+C.s[qc]+C.w[qc]+C.e[qc]+C.t[qc]+C.b[qc]);
        }
    }
}

//  ------------------------------------------------------------- line solver
//  Symmetric sweep; the tridiagonal solve in k is exact, so vertical
//  stretching costs nothing.  Identity rows sit in the system harmlessly
//  as diag 1, sub and super 0.
void reefmg_core::line_gs(sc_level &L,int l,int sweeps,int dir)
{
    const int nz=L.nz;
    std::vector<double> d(nz),rhs(nz),c(nz);

    const int p0=(dir==1?1:0), p1=(dir==0?1:2);

    for(int sw=0;sw<sweeps;++sw)
    for(int pass=p0;pass<p1;++pass)
    {
        halo(L);

        const int i0=(pass==0?0:L.nx-1);
        const int i1=(pass==0?L.nx:-1);
        const int di=(pass==0?1:-1);

        for(int i=i0;i!=i1;i+=di)
        for(int j=0;j<L.ny;++j)
        {
            const long col=L.idx(i,j,0);
            bool any=false;
            for(int k=0;k<nz;++k) if(L.act[col+k]){any=true; break;}
            if(!any) continue;

            for(int k=0;k<nz;++k)
            {
                const long q=col+k;
                double v=L.f[q];
                v-=L.n[q]*L.u[L.idx(i+1,j,k)];
                v-=L.s[q]*L.u[L.idx(i-1,j,k)];
                v-=L.w[q]*L.u[L.idx(i,j+1,k)];
                v-=L.e[q]*L.u[L.idx(i,j-1,k)];
                rhs[k]=v; d[k]=L.p[q];
            }

            c[0]=L.t[col]/d[0];
            rhs[0]/=d[0];
            for(int k=1;k<nz;++k)
            {
                const long q=col+k;
                const double m=d[k]-L.b[q]*c[k-1];
                c[k]=L.t[q]/m;
                rhs[k]=(rhs[k]-L.b[q]*rhs[k-1])/m;
            }
            L.u[col+nz-1]=rhs[nz-1];
            for(int k=nz-2;k>=0;--k)
            L.u[col+k]=rhs[k]-c[k]*L.u[col+k+1];
        }
    }
}

void reefmg_core::residual(sc_level &L)
{
    halo(L);

    for(int i=0;i<L.nx;++i)
    for(int j=0;j<L.ny;++j)
    for(int k=0;k<L.nz;++k)
    {
        const long q=L.idx(i,j,k);
        double v=L.p[q]*L.u[q]
                +L.n[q]*L.u[L.idx(i+1,j,k)]
                +L.s[q]*L.u[L.idx(i-1,j,k)]
                +L.w[q]*L.u[L.idx(i,j+1,k)]
                +L.e[q]*L.u[L.idx(i,j-1,k)];
        if(k<L.nz-1) v+=L.t[q]*L.u[q+1];
        if(k>0     ) v+=L.b[q]*L.u[q-1];
        L.r[q]=L.f[q]-v;
    }
}

void reefmg_core::restrict_xy(sc_level &F,sc_level &C)
{
    std::fill(C.f.begin(),C.f.end(),0.0);
    std::fill(C.u.begin(),C.u.end(),0.0);

    for(int I=0;I<C.nx;++I)
    for(int J=0;J<C.ny;++J)
    for(int k=0;k<C.nz;++k)
    {
        double s=0.0; int na=0;
        int nc=0;
        for(int a=0;a<C.rx;++a) for(int b=0;b<C.ry;++b)
        {
            if(C.rx*I+a>=F.nx || C.ry*J+b>=F.ny) continue;
            ++nc;
            const long qf=F.idx(C.rx*I+a,C.ry*J+b,k);
            if(F.act[qf]){s+=F.r[qf]; ++na;}
        }
        C.f[C.idx(I,J,k)] = (na>0 && nc>0)? s/double(nc) : 0.0;
    }
}

//  cell-centred bilinear in the plane, identity in z
void reefmg_core::prolong_xy(const sc_level &C,sc_level &F)
{
    for(int i=0;i<F.nx;++i)
    for(int j=0;j<F.ny;++j)
    {
        const int I=i/C.rx, J=j/C.ry;
        int Ia=I, Ja=J;
        double wx=1.0, wy=1.0;

        if(C.rx==2){Ia=(i%2==0)? I-1 : I+1; wx=0.75;}
        if(C.ry==2){Ja=(j%2==0)? J-1 : J+1; wy=0.75;}
        Ia=std::max(0,std::min(C.nx-1,Ia));
        Ja=std::max(0,std::min(C.ny-1,Ja));

        const double w00=wx*wy, w10=(1.0-wx)*wy, w01=wx*(1.0-wy), w11=(1.0-wx)*(1.0-wy);

        for(int k=0;k<F.nz;++k)
        {
            const long qf=F.idx(i,j,k);
            if(F.act[qf]==0) continue;
            F.u[qf] += w00*C.u[C.idx(I ,J ,k)] + w10*C.u[C.idx(Ia,J ,k)]
                     + w01*C.u[C.idx(I ,Ja,k)] + w11*C.u[C.idx(Ia,Ja,k)];
        }
    }
}

void reefmg_core::vcycle(int l,int pre,int post)
{
    sc_level &L=lev[l];

    if(l==(int)lev.size()-1)
    {
        line_gs(L,l,coarse_sweeps,2);
        return;
    }

    //  Alternating forward/backward costs one pass per smoothing step instead
    //  of two, and on a grid-aligned anisotropy it smooths just as well.
    const int dpre =(sweepstyle==1)?2:0;
    const int dpost=(sweepstyle==1)?2:1;

    line_gs(L,l,pre,dpre);
    residual(L);
    restrict_xy(L,lev[l+1]);
    vcycle(l+1,pre,post);
    prolong_xy(lev[l+1],L);
    line_gs(L,l,post,dpost);
}

double reefmg_core::dot(const sc_level &L,const std::vector<double> &a,
                                               const std::vector<double> &b) const
{
    double s=0.0;
    for(int i=0;i<L.nx;++i)
    for(int j=0;j<L.ny;++j)
    for(int k=0;k<L.nz;++k)
    {
        const long q=L.idx(i,j,k);
        if(L.act[q]) s+=a[q]*b[q];
    }
    double g=0.0;
    MPI_Allreduce(&s,&g,1,MPI_DOUBLE,MPI_SUM,comm);
    return g;
}

void reefmg_core::apply(sc_level &L,int l,const std::vector<double> &x,
                             std::vector<double> &y)
{
    std::vector<double> save;
    save.swap(L.u);
    L.u=x;
    halo(L);

    for(int i=0;i<L.nx;++i)
    for(int j=0;j<L.ny;++j)
    for(int k=0;k<L.nz;++k)
    {
        const long q=L.idx(i,j,k);
        double v=L.p[q]*L.u[q]
                +L.n[q]*L.u[L.idx(i+1,j,k)]
                +L.s[q]*L.u[L.idx(i-1,j,k)]
                +L.w[q]*L.u[L.idx(i,j+1,k)]
                +L.e[q]*L.u[L.idx(i,j-1,k)];
        if(k<L.nz-1) v+=L.t[q]*L.u[q+1];
        if(k>0     ) v+=L.b[q]*L.u[q-1];
        y[q]=v;
    }
    L.u.swap(save);
}

void reefmg_core::precondition(const std::vector<double> &rhs,
                                    std::vector<double> &x,int pre,int post)
{
    sc_level &F=lev[0];
    std::vector<double> savef,saveu;
    savef.swap(F.f); saveu.swap(F.u);
    F.f=rhs;
    F.u.assign(F.size(),0.0);
    vcycle(0,pre,post);
    x=F.u;
    F.f.swap(savef); F.u.swap(saveu);
}

int reefmg_core::solve_vcycle(double tol,int maxiter,double &relres,int pre,int post)
{
    sc_level &F=lev[0];
    const long N=F.size();
    for(long q=0;q<N;++q) if(F.act[q]==0){F.u[q]=0.0; F.f[q]=0.0;}

    const double bn=sqrt(dot(F,F.f,F.f));
    if(bn==0.0){relres=0.0; return 0;}

    int it=0;
    residual(F);
    double rn=sqrt(dot(F,F.r,F.r));

    while(rn/bn>tol && it<maxiter)
    {
        vcycle(0,pre,post);
        ++it;
        residual(F);
        rn=sqrt(dot(F,F.r,F.r));
    }
    relres=rn/bn;
    return it;
}

//  BiCGStab, V-cycle preconditioned.  The FNPF matrix is not symmetric, so a
//  Krylov wrapper is the safe default even though the cycle alone usually
//  converges.
int reefmg_core::solve(double tol,int maxiter,double &relres,int pre,int post)
{
    sc_level &F=lev[0];
    const long N=F.size();

    for(long q=0;q<N;++q) if(F.act[q]==0){F.u[q]=0.0; F.f[q]=0.0;}

    const double bn=sqrt(dot(F,F.f,F.f));
    if(bn==0.0){relres=0.0; std::fill(F.u.begin(),F.u.end(),0.0); return 0;}

    residual(F);
    kr=F.r;
    krhat=kr;

    double rho=1.0,alpha=1.0,omega=1.0;
    std::fill(kv.begin(),kv.end(),0.0);
    std::fill(kp.begin(),kp.end(),0.0);

    double rn=sqrt(dot(F,kr,kr));
    int it=0;

    while(rn/bn>tol && it<maxiter)
    {
        const double rho1=dot(F,krhat,kr);
        if(fabs(rho1)<1.0e-300) break;

        if(it==0)
        kp=kr;
        else
        {
            const double beta=(rho1/rho)*(alpha/omega);
            for(long q=0;q<N;++q) kp[q]=kr[q]+beta*(kp[q]-omega*kv[q]);
        }
        rho=rho1;

        precondition(kp,ky,pre,post);
        apply(F,0,ky,kv);

        const double den=dot(F,krhat,kv);
        if(fabs(den)<1.0e-300) break;
        alpha=rho/den;

        for(long q=0;q<N;++q) ks[q]=kr[q]-alpha*kv[q];

        //  standard half-step exit: if s is already small enough the second
        //  preconditioner application of this iteration is not needed
        const double sn=sqrt(dot(F,ks,ks));
        if(sn/bn<=tol)
        {
            for(long q=0;q<N;++q) F.u[q]+=alpha*ky[q];
            kr=ks;
            rn=sn;
            ++it;
            break;
        }

        precondition(ks,kz,pre,post);
        apply(F,0,kz,kt);

        const double tt=dot(F,kt,kt);
        omega=(tt>0.0)? dot(F,kt,ks)/tt : 0.0;

        for(long q=0;q<N;++q)
        {
            F.u[q]+=alpha*ky[q]+omega*kz[q];
            kr[q] =ks[q]-omega*kt[q];
        }

        ++it;
        rn=sqrt(dot(F,kr,kr));

        if(fabs(omega)<1.0e-300) break;
    }

    relres=rn/bn;
    return it;
}

//  V-cycles first, Krylov only if they misbehave.  A stalling cycle means the
//  coarse operator has lost too much of the fine one - steep bathymetry or a
//  strongly one-sided boundary row - and BiCGStab recovers it from wherever
//  the cycles left off.
int reefmg_core::solve_auto(double tol,int maxiter,double &relres,
                                 int pre,int post,int mode)
{
    if(mode==1)
    return solve(tol,maxiter,relres,pre,post);

    sc_level &F=lev[0];
    const long N=F.size();
    for(long q=0;q<N;++q) if(F.act[q]==0){F.u[q]=0.0; F.f[q]=0.0;}

    const double bn=sqrt(dot(F,F.f,F.f));
    if(bn==0.0){relres=0.0; std::fill(F.u.begin(),F.u.end(),0.0); return 0;}

    residual(F);
    double rn=sqrt(dot(F,F.r,F.r)), prev=rn;
    int it=0, bad=0;

    while(rn/bn>tol && it<maxiter)
    {
        vcycle(0,pre,post);
        ++it;
        residual(F);
        prev=rn;
        rn=sqrt(dot(F,F.r,F.r));

        if(rn>0.6*prev) ++bad; else bad=0;

        if(mode==2 && bad>=2)
        {
            ++nfallback;
            double rr;
            const int extra=solve(tol,maxiter-it,rr,pre,post);
            relres=rr;
            return it+extra;
        }
    }

    relres=rn/bn;
    return it;
}
