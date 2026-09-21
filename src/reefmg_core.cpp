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

namespace
{
    //  Everything the V-cycle kernels read, as raw pointers of type C.  The
    //  kernels are written once and instantiated for double and float; the
    //  arithmetic is always done in double, so float only changes how many
    //  bytes each coefficient costs to load.
    template<class C> struct sc_view
    {
        const C *p,*n,*s,*w,*e,*t,*b,*tc,*ti;
    };

    template<class C> sc_view<C> coef(const sc_level &L);

    template<> inline sc_view<double> coef<double>(const sc_level &L)
    {
        sc_view<double> v={L.p.data(),L.n.data(),L.s.data(),L.w.data(),L.e.data(),
                           L.t.data(),L.b.data(),L.tc.data(),L.ti.data()};
        return v;
    }

    template<> inline sc_view<float> coef<float>(const sc_level &L)
    {
        sc_view<float> v={L.pf.data(),L.nf.data(),L.sf.data(),L.wf.data(),L.ef.data(),
                          L.tf.data(),L.bf.data(),L.tcf.data(),L.tif.data()};
        return v;
    }
}

reefmg_core::reefmg_core()
{
    comm=MPI_COMM_NULL;
    myrank=nprocs=0;
    nbx0=nbx1=nby0=nby1=MPI_PROC_NULL;
    coarse_sweeps=16;
    nfallback=0;
    sweepstyle=0;
    pcbits=64;
    ordering=0;
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
void reefmg_core::halo_vec(sc_level &L,std::vector<double> &u)
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
    double *sx0=&u[L.idx(0,0,0)];
    double *sx1=&u[L.idx(L.nx-1,0,0)];
    double *rx0=&u[L.idx(-1,0,0)];
    double *rx1=&u[L.idx(L.nx,0,0)];

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
        p0[i*nz+k]=u[L.idx(i,0,k)];
        p1[i*nz+k]=u[L.idx(i,L.ny-1,k)];
    }

    MPI_Irecv(q0,cnty,MPI_DOUBLE,nby0,ty0,comm,&req[nreq++]);
    MPI_Irecv(q1,cnty,MPI_DOUBLE,nby1,ty1,comm,&req[nreq++]);
    MPI_Isend(p0,cnty,MPI_DOUBLE,nby0,ty1,comm,&req[nreq++]);
    MPI_Isend(p1,cnty,MPI_DOUBLE,nby1,ty0,comm,&req[nreq++]);
    MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);

    if(nby0!=MPI_PROC_NULL)
    for(int i=0;i<L.nx;++i) for(int k=0;k<nz;++k) u[L.idx(i,-1,k)]=q0[i*nz+k];
    if(nby1!=MPI_PROC_NULL)
    for(int i=0;i<L.nx;++i) for(int k=0;k<nz;++k) u[L.idx(i,L.ny,k)]=q1[i*nz+k];
}

void reefmg_core::halo(sc_level &L)
{
    halo_vec(L,L.u);
}

//  --------------------------------------------------------- coarse operators
//  Horizontal: a coarse face is two fine faces, so the coarse coefficient is
//  an eighth of their sum (half the conductance over four times the volume).
//  Vertical: a coarse face is four fine faces over four times the volume, so
//  the coefficient is their mean.  The row sum is carried across unchanged,
//  which is what keeps the free-surface Dirichlet term alive on every level.
void reefmg_core::coarsen()
{
    const bool fp32=(pcbits==32);

    if(fp32)
    for(size_t l=0;l<lev.size();++l)
    {
        sc_level &L=lev[l];
        const long N=L.size();
        if((long)L.pf.size()!=N)
        {
            L.pf.assign(N,0.0f); L.nf.assign(N,0.0f); L.sf.assign(N,0.0f);
            L.wf.assign(N,0.0f); L.ef.assign(N,0.0f); L.tf.assign(N,0.0f);
            L.bf.assign(N,0.0f); L.tcf.assign(N,0.0f); L.tif.assign(N,0.0f);
        }
    }

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

                //  Every fine cell is visited exactly once here, with all
                //  seven coefficients about to be loaded anyway, so the fp32
                //  mirror costs stores only - no separate pass over the
                //  hierarchy.  Identity rows mirror as identity rows.
                if(fp32)
                {
                    F.pf[qf]=(float)F.p[qf]; F.nf[qf]=(float)F.n[qf];
                    F.sf[qf]=(float)F.s[qf]; F.wf[qf]=(float)F.w[qf];
                    F.ef[qf]=(float)F.e[qf]; F.tf[qf]=(float)F.t[qf];
                    F.bf[qf]=(float)F.b[qf];
                }

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

    if(fp32)
    to_fp32();          // coarsest level only - the others were mirrored above

    factor_lines();
}

//  The column matrices are fixed for the whole solve - every sweep of every
//  V-cycle of every Krylov iteration sees the same ones - so their LU
//  factorisation is computed once here.  The smoother then needs only
//  multiplications; the two divisions per cell it used to perform on every
//  sweep were what the serial recurrence stalled on.
void reefmg_core::factor_lines()
{
    const bool fp32=(pcbits==32);

    for(size_t l=0;l<lev.size();++l)
    {
        sc_level &L=lev[l];
        const int nz=L.nz;
        const long N=L.size();

        if(!fp32 && (long)L.tc.size()!=N)
        {
            L.tc.assign(N,0.0);
            L.ti.assign(N,0.0);
        }
        L.colact.assign((long)L.nx*L.ny,0);
        L.zcol[0].clear();
        L.zcol[1].clear();

        for(int i=0;i<L.nx;++i)
        for(int j=0;j<L.ny;++j)
        {
            const long col=L.idx(i,j,0);

            char any=0;
            for(int k=0;k<nz;++k) if(L.act[col+k]){any=1; break;}
            L.colact[(long)i*L.ny+j]=any;
            if(any) L.zcol[(i+j)&1].push_back(col);

            //  The elimination runs in double whatever the storage precision;
            //  only the stored factors are rounded.  In fp32 mode the double
            //  factors are never needed, so they are neither kept nor allocated.
            double inv=1.0/L.p[col];
            double tcp=L.t[col]*inv;
            if(fp32){L.tif[col]=(float)inv; L.tcf[col]=(float)tcp;}
            else    {L.ti [col]=inv;        L.tc [col]=tcp;}

            for(int k=1;k<nz;++k)
            {
                const long q=col+k;
                inv=1.0/(L.p[q]-L.b[q]*tcp);
                tcp=L.t[q]*inv;
                if(fp32){L.tif[q]=(float)inv; L.tcf[q]=(float)tcp;}
                else    {L.ti [q]=inv;        L.tc [q]=tcp;}
            }
        }
    }
}

//  ------------------------------------------------------------- line solver
//  Symmetric sweep; the tridiagonal solve in k is exact, so vertical
//  stretching costs nothing.  Identity rows sit in the system harmlessly
//  as diag 1, sub and super 0.
void reefmg_core::line_gs(sc_level &L,int l,int sweeps,int dir)
{
    if(ordering==1)
    {
        if(pcbits==32) line_zebra_t<float >(L,sweeps,dir);
        else           line_zebra_t<double>(L,sweeps,dir);
    }
    else
    {
        if(pcbits==32) line_gs_t<float >(L,sweeps,dir);
        else           line_gs_t<double>(L,sweeps,dir);
    }
}

//  Red-black line Gauss-Seidel with the column solves batched across SIMD
//  lanes.  Every red column's horizontal neighbours are black and vice versa,
//  so all columns of one colour are independent and can be solved together.
//  ZB columns are gathered into a transposed scratch block [k][lane]; the
//  tridiagonal recurrence then runs down k with all lanes in step, which is
//  what the compiler vectorises.  The block holds five arrays of nz*ZB
//  doubles and stays in L1 for any realistic number of layers.
//
//  Pre-smoothing solves red then black, post-smoothing black then red, so the
//  V-cycle keeps the same forward/backward symmetry as the lexicographic
//  sweep.  A halo exchange precedes each half-sweep so each colour sees the
//  other colour's latest values across process boundaries; colours are local,
//  so a pair of same-coloured columns across a boundary couples block-Jacobi
//  style, exactly as every interface did with the lexicographic sweep.
static const int ZB=8;

template<class C>
void reefmg_core::line_zebra_t(sc_level &L,int sweeps,int dir)
{
    const sc_view<C> V=coef<C>(L);
    const int nz=L.nz;
    const long sx=(long)(L.ny+2)*nz, sy=nz;

    const long ns=(long)nz*ZB;
    if((long)zr.size()<ns){zr.resize(ns); zb.resize(ns); zti.resize(ns); ztc.resize(ns);}
    double *R=&zr[0], *Bb=&zb[0], *Ti=&zti[0], *Tc=&ztc[0];

    const double *f=&L.f[0];
    double *u=&L.u[0];

    for(int sw=0;sw<sweeps;++sw)
    for(int h=0;h<2;++h)
    {
        const int colour=(dir==1)? 1-h : h;
        const std::vector<long> &cols=L.zcol[colour];
        const long nc=cols.size();

        halo(L);

        for(long c0=0;c0<nc;c0+=ZB)
        {
            //  a short final batch repeats its last column; the duplicates
            //  read identical inputs and write identical values
            long cb[ZB];
            for(int b=0;b<ZB;++b) cb[b]=cols[std::min(c0+b,nc-1)];

            //  gather: right-hand side and factors, transposed to [k][lane]
            for(int b=0;b<ZB;++b)
            {
                const long c=cb[b];
                for(int k=0;k<nz;++k)
                {
                    const long q=c+k;
                    R [k*ZB+b]=f[q]-V.n[q]*u[q+sx]-V.s[q]*u[q-sx]
                                   -V.w[q]*u[q+sy]-V.e[q]*u[q-sy];
                    Bb[k*ZB+b]=V.b[q];
                    Ti[k*ZB+b]=V.ti[q];
                    Tc[k*ZB+b]=V.tc[q];
                }
            }

            //  forward and back substitution, all lanes in step
            for(int b=0;b<ZB;++b) R[b]*=Ti[b];
            for(int k=1;k<nz;++k)
            {
                double *r=R+k*ZB; const double *rm=R+(k-1)*ZB;
                const double *bb=Bb+k*ZB, *ti=Ti+k*ZB;
                for(int b=0;b<ZB;++b) r[b]=(r[b]-bb[b]*rm[b])*ti[b];
            }
            for(int k=nz-2;k>=0;--k)
            {
                double *r=R+k*ZB; const double *rp=R+(k+1)*ZB;
                const double *tc=Tc+k*ZB;
                for(int b=0;b<ZB;++b) r[b]-=tc[b]*rp[b];
            }

            //  scatter
            for(int b=0;b<ZB;++b)
            {
                double *uc=u+cb[b];
                for(int k=0;k<nz;++k) uc[k]=R[k*ZB+b];
            }
        }
    }
}

template<class C>
void reefmg_core::line_gs_t(sc_level &L,int sweeps,int dir)
{
    const sc_view<C> V=coef<C>(L);
    const int nz=L.nz;
    std::vector<double> rhs(nz);
    double *r=&rhs[0];

    const double *Lf=&L.f[0];
    const C *Ln=V.n, *Ls=V.s, *Lw=V.w, *Le=V.e, *Lb=V.b, *Tc=V.tc, *Ti=V.ti;
    double *Lu=&L.u[0];

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
            if(L.colact[(long)i*L.ny+j]==0) continue;

            const long col=L.idx(i,j,0);
            const double *un=Lu+L.idx(i+1,j,0), *us=Lu+L.idx(i-1,j,0);
            const double *uw=Lu+L.idx(i,j+1,0), *ue=Lu+L.idx(i,j-1,0);

            //  right-hand side: contiguous in k
            for(int k=0;k<nz;++k)
            {
                const long q=col+k;
                r[k]=Lf[q]-Ln[q]*un[k]-Ls[q]*us[k]-Lw[q]*uw[k]-Le[q]*ue[k];
            }

            //  forward and back substitution with the stored factorisation
            r[0]*=Ti[col];
            for(int k=1;k<nz;++k)
            r[k]=(r[k]-Lb[col+k]*r[k-1])*Ti[col+k];

            double *uc=Lu+col;
            uc[nz-1]=r[nz-1];
            for(int k=nz-2;k>=0;--k)
            uc[k]=r[k]-Tc[col+k]*uc[k+1];
        }
    }
}


//  y = A*u on one column, split so every loop is contiguous in k and free
//  of branches; the two vertical couplings are applied in their own loops.
template<class C>
static inline void column_apply(const sc_view<C> &V,int nz,const double *u,long col,
                                long cn,long cs,long cw,long ce,double *y)
{
    const C *P=V.p+col, *N=V.n+col, *S=V.s+col, *W=V.w+col, *E=V.e+col;
    const C *T=V.t+col, *B=V.b+col;
    const double *uc=u+col, *un=u+cn, *us=u+cs, *uw=u+cw, *ue=u+ce;

    for(int k=0;k<nz;++k)
    y[k]=P[k]*uc[k]+N[k]*un[k]+S[k]*us[k]+W[k]*uw[k]+E[k]*ue[k];

    for(int k=0;k<nz-1;++k) y[k]+=T[k]*uc[k+1];
    for(int k=1;k<nz;  ++k) y[k]+=B[k]*uc[k-1];
}

void reefmg_core::residual(sc_level &L)
{
    residual_t<double>(L);
}

void reefmg_core::residual_cycle(sc_level &L)
{
    if(pcbits==32) residual_t<float >(L);
    else           residual_t<double>(L);
}

template<class C>
void reefmg_core::residual_t(sc_level &L)
{
    const sc_view<C> V=coef<C>(L);
    halo(L);

    for(int i=0;i<L.nx;++i)
    for(int j=0;j<L.ny;++j)
    {
        const long col=L.idx(i,j,0);
        double *r=&L.r[col];
        column_apply<C>(V,L.nz,&L.u[0],col,L.idx(i+1,j,0),L.idx(i-1,j,0),
                                          L.idx(i,j+1,0),L.idx(i,j-1,0),r);
        const double *f=&L.f[col];
        for(int k=0;k<L.nz;++k) r[k]=f[k]-r[k];
    }
}

//  Mirror the coarsest level into single precision.  Every other level is
//  mirrored during coarsen(), while its coefficients are loaded to build the
//  next level; the coarsest is never a fine level, so it is done here.  It is
//  also the only level when the hierarchy has not been coarsened at all.
void reefmg_core::to_fp32()
{
    sc_level &L=lev.back();
    const long N=L.size();

    for(long q=0;q<N;++q)
    {
        L.pf[q]=(float)L.p[q]; L.nf[q]=(float)L.n[q]; L.sf[q]=(float)L.s[q];
        L.wf[q]=(float)L.w[q]; L.ef[q]=(float)L.e[q]; L.tf[q]=(float)L.t[q];
        L.bf[q]=(float)L.b[q];
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
    residual_cycle(L);
    restrict_xy(L,lev[l+1]);
    vcycle(l+1,pre,post);
    prolong_xy(lev[l+1],L);
    line_gs(L,l,post,dpost);
}

//  Inactive entries are zero in every vector this is called with - identity
//  rows carry f=0 and u=0, and the Krylov updates preserve that - so the
//  activity test the loop used to make per cell is unnecessary.  Four partial
//  sums break the single dependency chain of one running total.
double reefmg_core::dot(const sc_level &L,const std::vector<double> &a,
                                               const std::vector<double> &b) const
{
    const int nz=L.nz;
    double s0=0.0,s1=0.0,s2=0.0,s3=0.0;

    for(int i=0;i<L.nx;++i)
    for(int j=0;j<L.ny;++j)
    {
        const double *pa=&a[L.idx(i,j,0)], *pb=&b[L.idx(i,j,0)];
        int k=0;
        for(;k+3<nz;k+=4)
        {
            s0+=pa[k  ]*pb[k  ];
            s1+=pa[k+1]*pb[k+1];
            s2+=pa[k+2]*pb[k+2];
            s3+=pa[k+3]*pb[k+3];
        }
        for(;k<nz;++k) s0+=pa[k]*pb[k];
    }

    double s=(s0+s1)+(s2+s3), g=0.0;
    MPI_Allreduce(&s,&g,1,MPI_DOUBLE,MPI_SUM,comm);
    return g;
}

//  x is exchanged in place.  The previous version swapped L.u out and
//  assigned x into it, which allocated, copied and freed a full-size vector
//  on every product - and for large vectors the allocation goes through mmap,
//  so every call also page-faulted the buffer in from scratch.
void reefmg_core::apply(sc_level &L,int l,std::vector<double> &x,
                             std::vector<double> &y)
{
    const sc_view<double> V64=coef<double>(L);
    halo_vec(L,x);

    for(int i=0;i<L.nx;++i)
    for(int j=0;j<L.ny;++j)
    {
        const long col=L.idx(i,j,0);
        column_apply<double>(V64,L.nz,&x[0],col,L.idx(i+1,j,0),L.idx(i-1,j,0),
                                                L.idx(i,j+1,0),L.idx(i,j-1,0),&y[col]);
    }
}

//  Buffers are swapped rather than copied: rhs becomes the fine-level
//  right-hand side and x the fine-level solution for the duration of the
//  cycle, then both are handed back.  The cycle does not modify the fine f.
void reefmg_core::precondition(std::vector<double> &rhs,
                                    std::vector<double> &x,int pre,int post)
{
    sc_level &F=lev[0];
    F.f.swap(rhs);
    F.u.swap(x);
    std::fill(F.u.begin(),F.u.end(),0.0);
    vcycle(0,pre,post);
    F.u.swap(x);
    F.f.swap(rhs);
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
