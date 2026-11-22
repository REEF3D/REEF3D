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

#include"nhflow_si_nhmg.h"
#include<cmath>
#include<cstring>
#include<algorithm>
#include<iostream>
#include<cstdlib>

namespace
{
const int dix[5] = {0, 1,-1, 0, 0};
const int diy[5] = {0, 0, 0, 1,-1};

template<class T> MPI_Datatype mpitype();
template<> MPI_Datatype mpitype<float>() { return MPI_FLOAT; }
template<> MPI_Datatype mpitype<double>() { return MPI_DOUBLE; }

inline int hix(int ny, int ix, int iy) { return (ix+1)*(ny+2) + (iy+1); }

// the halo of one column of a vector with halo (nx+2) x (ny+2) x B from the ranks nb (-x, +x, -y, +y); one rank
// without neighbours: the halo stays zero (as with MPI_PROC_NULL)
template<class T> void halo(MPI_Comm comm, const int nb[4], bool alone, int B, int nx, int ny, std::vector<T> &v,
                            std::vector<T> &sbuf, std::vector<T> &rbuf)
{
    if(alone)
    return;

    auto strip = [&](bool xdir, int idx, bool get, std::vector<T> &buf)
    {
        const int n = xdir ? ny : nx;
        buf.resize((size_t)n*B);

        for(int m=0; m<n; ++m)
        {
            const int ix = xdir ? idx : m, iy = xdir ? m : idx;
            T *p = &v[(size_t)hix(ny,ix,iy)*B];

            if(get)
            memcpy(&buf[(size_t)m*B],p,B*sizeof(T));
            else
            memcpy(p,&buf[(size_t)m*B],B*sizeof(T));
        }
    };

    auto pass = [&](bool xdir, int send_idx, int dest, int recv_idx, int src, int tag)
    {
        const int n = (xdir ? ny : nx)*B;
        strip(xdir,send_idx,true,sbuf);
        rbuf.assign(n,T(0));
        MPI_Sendrecv(sbuf.data(),n,mpitype<T>(),dest,tag,rbuf.data(),n,mpitype<T>(),src,tag,comm,MPI_STATUS_IGNORE);

        if(src!=MPI_PROC_NULL)
        strip(xdir,recv_idx,false,rbuf);
    };

    pass(true,0,nb[0],nx,nb[1],7610);
    pass(true,nx-1,nb[1],-1,nb[0],7611);

    if(nb[2]!=MPI_PROC_NULL || nb[3]!=MPI_PROC_NULL || ny>1)
    {
    pass(false,0,nb[2],ny,nb[3],7612);
    pass(false,ny-1,nb[3],-1,nb[2],7613);
    }
}

// y[r] += sign sum_c a_rc x[c] over the nonzero pattern pat[r] of a block (a_rc stored at a[o])
template<class T, class U> inline void block_mul(int B, const std::vector<nhmg_pent> *pat, const T *a, const U *x, U *y, U sign)
{
    for(int r=0; r<B; ++r)
    {
        U t = 0;
        for(const nhmg_pent &e : pat[r])
        t += (U)a[e.o]*x[e.c];
        y[r] += sign*t;
    }
}


// the blocks of the coupled system with the mask of nh_mask (layout "arrow"): the diagonal block in full (row-major),
// a neighbour block row 0 in full (B entries) and the rows r >= 1 at c = 0, r, r+1 (3 entries, the last row 2);
// fixed-size kernels without the pattern lists
template<int NB, class T, class U> inline void arrow_nb(const T *a, const U *x, U *y, U sign)
{
    U t = 0;
    for(int c=0; c<NB; ++c)
    t += (U)a[c]*x[c];
    y[0] += sign*t;
    const T *e = a + NB;
    for(int r=1; r<NB-1; ++r, e+=3)
    y[r] += sign*((U)e[0]*x[0] + (U)e[1]*x[r] + (U)e[2]*x[r+1]);
    if(NB>1)
    y[NB-1] += sign*((U)e[0]*x[0] + (U)e[1]*x[NB-1]);
}
template<int NB, class T, class U> inline void dense_b(const T *a, const U *x, U *y, U sign)
{
    for(int r=0; r<NB; ++r)
    {
        U t = 0;
        for(int c=0; c<NB; ++c)
        t += (U)a[r*NB+c]*x[c];
        y[r] += sign*t;
    }
}
template<class T, class U> inline void arrow_nb_rt(int B, const T *a, const U *x, U *y, U sign)
{
    U t = 0;
    for(int c=0; c<B; ++c)
    t += (U)a[c]*x[c];
    y[0] += sign*t;
    const T *e = a + B;
    for(int r=1; r<B-1; ++r, e+=3)
    y[r] += sign*((U)e[0]*x[0] + (U)e[1]*x[r] + (U)e[2]*x[r+1]);
    if(B>1)
    y[B-1] += sign*((U)e[0]*x[0] + (U)e[1]*x[B-1]);
}
template<class T, class U> inline void dense_rt(int B, const T *a, const U *x, U *y, U sign)
{
    for(int r=0; r<B; ++r)
    {
        U t = 0;
        for(int c=0; c<B; ++c)
        t += (U)a[r*B+c]*x[c];
        y[r] += sign*t;
    }
}
// block d of a column: arrow layout by the fixed kernels, otherwise by the pattern
template<class T, class U> inline void bmul(int B, bool arrow, int d, const std::vector<nhmg_pent> *pat, const T *a, const U *x, U *y, U sign)
{
    if(!arrow)
    { block_mul<T,U>(B,pat,a,x,y,sign); return; }

    if(d==0)
    switch(B)
    {
    case 2: dense_b<2,T,U>(a,x,y,sign); return;
    case 3: dense_b<3,T,U>(a,x,y,sign); return;
    case 4: dense_b<4,T,U>(a,x,y,sign); return;
    case 5: dense_b<5,T,U>(a,x,y,sign); return;
    case 6: dense_b<6,T,U>(a,x,y,sign); return;
    case 7: dense_b<7,T,U>(a,x,y,sign); return;
    case 8: dense_b<8,T,U>(a,x,y,sign); return;
    case 9: dense_b<9,T,U>(a,x,y,sign); return;
    case 11: dense_b<11,T,U>(a,x,y,sign); return;
    default: dense_rt<T,U>(B,a,x,y,sign); return;
    }

    switch(B)
    {
    case 2: arrow_nb<2,T,U>(a,x,y,sign); return;
    case 3: arrow_nb<3,T,U>(a,x,y,sign); return;
    case 4: arrow_nb<4,T,U>(a,x,y,sign); return;
    case 5: arrow_nb<5,T,U>(a,x,y,sign); return;
    case 6: arrow_nb<6,T,U>(a,x,y,sign); return;
    case 7: arrow_nb<7,T,U>(a,x,y,sign); return;
    case 8: arrow_nb<8,T,U>(a,x,y,sign); return;
    case 9: arrow_nb<9,T,U>(a,x,y,sign); return;
    case 11: arrow_nb<11,T,U>(a,x,y,sign); return;
    default: arrow_nb_rt<T,U>(B,a,x,y,sign); return;
    }
}

// bmul with the block size NB fixed at compile time (NB = 0: B at run time), for the loops of dispatchB
template<int NB, class T, class U> inline void bmulB(int B, bool arrow, int d, const std::vector<nhmg_pent> *pat, const T *a, const U *x, U *y, U sign)
{
    if constexpr(NB==0)
    bmul<T,U>(B,arrow,d,pat,a,x,y,sign);

    else
    {
        if(!arrow)
        block_mul<T,U>(NB,pat,a,x,y,sign);
        else if(d==0)
        dense_b<NB,T,U>(a,x,y,sign);
        else
        arrow_nb<NB,T,U>(a,x,y,sign);
    }
}

// f.template operator()<NB>() with NB = B for the common column sizes (K = 1 .. 10 layers), otherwise NB = 0
template<class F> inline void dispatchB(int B, F &&f)
{
    switch(B)
    {
    case 2: f.template operator()<2>(); return;
    case 3: f.template operator()<3>(); return;
    case 4: f.template operator()<4>(); return;
    case 5: f.template operator()<5>(); return;
    case 6: f.template operator()<6>(); return;
    case 7: f.template operator()<7>(); return;
    case 8: f.template operator()<8>(); return;
    case 9: f.template operator()<9>(); return;
    case 11: f.template operator()<11>(); return;
    default: f.template operator()<0>(); return;
    }
}

// the inverse of a B x B block (LU with partial pivoting; a zero pivot: 1, as the LU of the column solves before);
// a: B*B + B scratch, pv: B
void invert(int B, const double *a0, double *inv, double *a, int *pv)
{
    std::copy(a0,a0+B*B,a);

    for(int k=0; k<B; ++k)
    {
        int pr = k;
        double pm = fabs(a[k*B+k]);

        for(int r=k+1; r<B; ++r)
        if(fabs(a[r*B+k])>pm) { pm = fabs(a[r*B+k]); pr = r; }

        pv[k] = pr;

        if(pr!=k)
        for(int c=0; c<B; ++c)
        std::swap(a[k*B+c],a[pr*B+c]);

        const double d = a[k*B+k];

        if(d==0.0)
        {
        a[k*B+k] = 1.0;
        continue;
        }

        for(int r=k+1; r<B; ++r)
        {
            const double f = a[r*B+k]/d;
            a[r*B+k] = f;

            if(f!=0.0)
            for(int c=k+1; c<B; ++c)
            a[r*B+c] -= f*a[k*B+c];
        }
    }

    double *x = a + B*B;
    for(int e=0; e<B; ++e)
    {
        std::fill(x,x+B,0.0);
        x[e] = 1.0;

        for(int k=0; k<B; ++k)
        if(pv[k]!=k)
        std::swap(x[k],x[pv[k]]);

        for(int k=0; k<B; ++k)
        for(int r=k+1; r<B; ++r)
        x[r] -= a[r*B+k]*x[k];

        for(int k=B-1; k>=0; --k)
        {
            double s = x[k];
            for(int c=k+1; c<B; ++c)
            s -= a[k*B+c]*x[c];
            x[k] = s/a[k*B+k];
        }

        for(int r=0; r<B; ++r)
        inv[r*B+e] = x[r];
    }
}
}

// ---- the preconditioner in the precision T
template<class T> class nhflow_si_nhmg_core
{
public:
    struct level
    {
        int nx=0, ny=0;
        std::vector<T> A;       // 5 blocks per column (CS entries, the storage of level 0)
        std::vector<T> Di;      // the inverse of the diagonal block
        std::vector<T> x, b, r;     // with the halo of one column
        std::vector<T> kb, v1, w1, v2, w2;     // K-cycle
        std::vector<std::vector<nhmg_pent>> pat;      // the columns c of row r of block d that can be nonzero: pat[d*B + r]
    };

    MPI_Comm comm;
    int B=0, BB=0, CS=0;
    int offs[5];
    int nb[4];
    bool alone = true;
    std::vector<level> L;
    std::vector<T> sbuf, rbuf, tmp, tmp2;
    const nhflow_si_nhmg *o = nullptr;     // the options

    void setup(MPI_Comm cm, int b, int nx, int ny, const int n[4], bool aln, const nhflow_si_nhmg *par)
    {
        comm = cm; B = b; BB = B*B; alone = aln; o = par;
        CS = o->CS;
        for(int q=0; q<5; ++q) offs[q] = o->offs[q];
        for(int q=0; q<4; ++q) nb[q] = n[q];

        int nl = 1, ax = nx, ay = ny;
        while((ax>2 || ay>2) && nl<30)
        {
        ax = (ax+1)/2;
        ay = (ay+1)/2;
        ++nl;
        }

        int gnl = nl;
        MPI_Allreduce(&nl,&gnl,1,MPI_INT,MPI_MIN,comm);

        L.assign(gnl,level());
        ax = nx; ay = ny;

        for(int l=0; l<gnl; ++l)
        {
        level &v = L[l];
        v.nx = ax;
        v.ny = ay;
        v.A.assign((size_t)ax*ay*CS,T(0));
        v.Di.assign((size_t)ax*ay*BB,T(0));
        v.x.assign((size_t)(ax+2)*(ay+2)*B,T(0));
        v.b = v.x;
        v.r = v.x;

        // the K-cycle vectors only on the levels it corrects (1 .. nl-2)
        if(l>0 && l<gnl-1)
        v.kb = v.v1 = v.w1 = v.v2 = v.w2 = v.x;
        ax = (ax+1)/2;
        ay = (ay+1)/2;
        }

        tmp.assign(B,T(0));
        tmp2.assign(B,T(0));
    }

    // the levels from the level-0 blocks (double) and their pattern: the Galerkin sums in double, the blocks and the
    // inverses of the diagonal blocks stored in T
    void build(const std::vector<double> &A0, const std::vector<std::vector<nhmg_pent>> &pat0)
    {
        // fine: level 0 (A0) and then the double sums of the coarse levels (fd, cd)
        std::vector<double> fd, cdv, inv(BB), sc(BB+B);
        std::vector<int> pv(B);
        const double *fine = A0.data();
        L[0].pat = pat0;

        for(int l=0; l<(int)L.size(); ++l)
        {
            level &f = L[l];
            const size_t nc = (size_t)f.nx*f.ny;

            for(size_t q=0; q<nc*CS; ++q)
            f.A[q] = (T)fine[q];

            // the diagonal block: in full at the start of the column
            for(size_t q=0; q<nc; ++q)
            {
            invert(B,&fine[q*CS],inv.data(),sc.data(),pv.data());
            for(int e=0; e<BB; ++e)
            f.Di[q*BB+e] = (T)inv[e];
            }

            if(l+1==(int)L.size())
            break;

            // the coarse operator: the blocks of the fine columns of an aggregate summed, a coupling within the aggregate
            // into its diagonal block; its pattern: block d (> 0) of the fine one, the diagonal block the union
            level &c = L[l+1];
            std::vector<double> &coarse = cdv;
            coarse.assign((size_t)c.nx*c.ny*CS,0.0);

            for(int ix=0; ix<f.nx; ++ix)
            for(int iy=0; iy<f.ny; ++iy)
            {
                const int Ix = ix/2, Iy = iy/2;

                for(int d=0; d<5; ++d)
                {
                    const double *a = &fine[(size_t)(ix*f.ny + iy)*CS + offs[d]];
                    const int jx = ix + dix[d], jy = iy + diy[d];
                    int dc;

                    if(jx<0 || jx>=f.nx || jy<0 || jy>=f.ny)
                    dc = d;

                    else
                    {
                    const int Jx = jx/2, Jy = jy/2;
                    dc = 0;
                    if(Jx==Ix+1) dc = 1;
                    if(Jx==Ix-1) dc = 2;
                    if(Jy==Iy+1) dc = 3;
                    if(Jy==Iy-1) dc = 4;
                    }

                    // the stored entries (r B + c, place) of block d; into the full diagonal block or the same place of
                    // block d (the same pattern)
                    double *t = &coarse[(size_t)(Ix*c.ny + Iy)*CS + offs[dc]];
                    if(dc==0)
                    for(const auto &e : o->sl[d])
                    t[e.first] += a[e.second];
                    else
                    for(const auto &e : o->sl[d])
                    t[e.second] += a[e.second];
                }
            }

            c.pat = f.pat;
            for(int r=0; r<B; ++r)
            {
                std::vector<char> u(B,0);
                for(int d=0; d<5; ++d)
                for(const nhmg_pent &e : f.pat[(size_t)d*B + r])
                u[e.c] = 1;

                std::vector<nhmg_pent> &p0 = c.pat[r];
                p0.clear();
                for(int cc=0; cc<B; ++cc)
                if(u[cc]) p0.push_back({cc,r*B+cc});
            }

            fd.swap(cdv);
            fine = fd.data();
        }
    }

    void exchange(level &l, std::vector<T> &v) { halo<T>(comm,nb,alone,B,l.nx,l.ny,v,sbuf,rbuf); }

    void smooth(int lev, int sweeps, bool zero=false)
    {
        // zero: x = 0 on entry, the first half-sweep has no neighbour products (the same result)
        level &v = L[lev];
        const bool ar = o->arrow;
        T *tp = tmp.data();

        for(int s=0; s<sweeps; ++s)
        for(int colour=0; colour<2; ++colour)
        {
            const bool first = zero && s==0 && colour==0;

            if(!first)
            exchange(v,v.x);

            // the columns of this colour only: (ix + iy) & 1 == colour
            dispatchB(B,[&]<int NB>()
            {
            for(int ix=0; ix<v.nx; ++ix)
            for(int iy=((ix^colour)&1); iy<v.ny; iy+=2)
            {
                const size_t q = (size_t)(ix*v.ny + iy);
                const T *bb = &v.b[(size_t)hix(v.ny,ix,iy)*B];

                for(int r=0; r<B; ++r)
                tp[r] = bb[r];

                if(!first)
                for(int d=1; d<5; ++d)
                {
                    const T *a = &v.A[q*CS + offs[d]];
                    const T *xx = &v.x[(size_t)hix(v.ny,ix+dix[d],iy+diy[d])*B];
                    bmulB<NB,T,T>(B,ar,d,&v.pat[(size_t)d*B],a,xx,tp,T(-1));
                }

                T *xx = &v.x[(size_t)hix(v.ny,ix,iy)*B];
                for(int r=0; r<B; ++r)
                xx[r] = 0;
                bmulB<NB,T,T>(B,true,0,nullptr,&v.Di[q*BB],tp,xx,T(1));
            }
            });
        }
    }

    // r = b - A x; skip >= 0: the colour relaxed last (its columns solved with the current neighbours: r = 0 there), the
    // columns at a rank face computed in full (the colours of the ranks need not match there)
    void residual(int lev, int skip=-1)
    {
        level &v = L[lev];
        const bool ar = o->arrow;
        exchange(v,v.x);

        dispatchB(B,[&]<int NB>()
        {
        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            const size_t q = (size_t)(ix*v.ny + iy);
            T *rr = &v.r[(size_t)hix(v.ny,ix,iy)*B];

            if(skip>=0 && ((ix+iy)&1)==skip && (alone || (ix>0 && ix<v.nx-1 && iy>0 && iy<v.ny-1)))
            {
            for(int r=0; r<B; ++r)
            rr[r] = T(0);
            continue;
            }

            const T *bb = &v.b[(size_t)hix(v.ny,ix,iy)*B];

            for(int r=0; r<B; ++r)
            rr[r] = bb[r];

            for(int d=0; d<5; ++d)
            {
                const T *a = &v.A[q*CS + offs[d]];
                const T *xx = &v.x[(size_t)hix(v.ny,ix+dix[d],iy+diy[d])*B];
                bmulB<NB,T,T>(B,ar,d,&v.pat[(size_t)d*B],a,xx,rr,T(-1));
            }
        }
        });
    }

    void cycle(int lev)
    {
        level &v = L[lev];
        std::fill(v.x.begin(),v.x.end(),T(0));

        if(lev==(int)L.size()-1)
        {
        smooth(lev,o->ncoarse,true);
        return;
        }

        smooth(lev,o->npre,true);
        residual(lev,(o->fused && o->npre>0) ? 1 : -1);

        level &c = L[lev+1];
        std::fill(c.b.begin(),c.b.end(),T(0));

        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            T *cb = &c.b[(size_t)hix(c.ny,ix/2,iy/2)*B];
            const T *rr = &v.r[(size_t)hix(v.ny,ix,iy)*B];
            for(int r=0; r<B; ++r)
            cb[r] += rr[r];
        }

        if(o->kcycle && lev+1<(int)L.size()-1)
        kcorrect(lev+1);
        else
        cycle(lev+1);

        const T w = (T)o->omega_c;
        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            T *xx = &v.x[(size_t)hix(v.ny,ix,iy)*B];
            const T *cx = &c.x[(size_t)hix(c.ny,ix/2,iy/2)*B];
            for(int r=0; r<B; ++r)
            xx[r] += w*cx[r];
        }

        smooth(lev,o->npost);
    }

    void apply(int lev, std::vector<T> &x, std::vector<T> &y)
    {
        level &v = L[lev];
        const bool ar = o->arrow;
        exchange(v,x);

        dispatchB(B,[&]<int NB>()
        {
        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            const size_t q = (size_t)(ix*v.ny + iy);
            T *yy = &y[(size_t)hix(v.ny,ix,iy)*B];
            std::fill(yy,yy+B,T(0));

            for(int d=0; d<5; ++d)
            {
                const T *a = &v.A[q*CS + offs[d]];
                const T *xx = &x[(size_t)hix(v.ny,ix+dix[d],iy+diy[d])*B];
                bmulB<NB,T,T>(B,ar,d,&v.pat[(size_t)d*B],a,xx,yy,T(1));
            }
        }
        });
    }

    void dots(int lev, const std::vector<T> *const *u, const std::vector<T> *const *w, int n, double *res)
    {
        level &v = L[lev];

        for(int k=0; k<n; ++k)
        {
            double s = 0.0;
            for(int ix=0; ix<v.nx; ++ix)
            for(int iy=0; iy<v.ny; ++iy)
            {
                const size_t h = (size_t)hix(v.ny,ix,iy)*B;
                for(int r=0; r<B; ++r)
                s += (double)(*u[k])[h+r]*(double)(*w[k])[h+r];
            }
            res[k] = s;
        }

        MPI_Allreduce(MPI_IN_PLACE,res,n,MPI_DOUBLE,MPI_SUM,comm);
    }

    // K-cycle (Notay): the coarse-level correction by two steps of GCR on A_lev e = b_lev, preconditioned by the cycle
    // of the level; the second step only when the first leaves more than ktol of |b|^2
    void kcorrect(int lev)
    {
        level &v = L[lev];
        v.kb = v.b;

        cycle(lev);
        v.v1 = v.x;
        apply(lev,v.v1,v.w1);

        double d[3];
        {
        const std::vector<T> *u[3] = {&v.w1,&v.w1,&v.kb};
        const std::vector<T> *w[3] = {&v.w1,&v.kb,&v.kb};
        dots(lev,u,w,3,d);
        }

        const double rho1 = d[0];
        if(rho1<=0.0)
        return;

        const double a1 = d[1]/rho1;
        const double rn2 = d[2] - d[1]*d[1]/rho1;     // |b - a1 w1|^2

        if(rn2 <= o->ktol*d[2])
        {
            for(size_t q=0; q<v.x.size(); ++q)
            v.x[q] = (T)(a1*v.v1[q]);
            return;
        }

        for(size_t q=0; q<v.b.size(); ++q)
        v.b[q] = (T)(v.kb[q] - a1*v.w1[q]);

        cycle(lev);
        v.v2 = v.x;
        apply(lev,v.v2,v.w2);

        {
        const std::vector<T> *u[3] = {&v.w1,&v.w2,&v.w2};
        const std::vector<T> *w[3] = {&v.w2,&v.w2,&v.b};
        dots(lev,u,w,3,d);
        }

        const double g = d[0]/rho1;
        const double rho2 = d[1] - d[0]*d[0]/rho1;            // |w2 - g w1|^2
        const double a2 = rho2>1.0e-300*d[1] ? d[2]/rho2 : 0.0;    // (w2 - g w1).r2 = w2.r2 (r2 orthogonal to w1)

        for(size_t q=0; q<v.x.size(); ++q)
        v.x[q] = (T)(a1*v.v1[q] + a2*(v.v2[q] - g*v.v1[q]));
    }

    void vcycle(const double *r, double *z)
    {
        level &v = L[0];
        std::fill(v.b.begin(),v.b.end(),T(0));

        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
        T *bb = &v.b[(size_t)hix(v.ny,ix,iy)*B];
        const double *rr = &r[(size_t)(ix*v.ny + iy)*B];
        for(int t=0; t<B; ++t)
        bb[t] = (T)rr[t];
        }

        cycle(0);

        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
        const T *xx = &v.x[(size_t)hix(v.ny,ix,iy)*B];
        double *zz = &z[(size_t)(ix*v.ny + iy)*B];
        for(int t=0; t<B; ++t)
        zz[t] = (double)xx[t];
        }
    }
};

template class nhflow_si_nhmg_core<float>;
template class nhflow_si_nhmg_core<double>;

// ---- level 0 and the interface
void nhflow_si_nhmg::setup(MPI_Comm cm, int b, int nx, int ny, const int n[4], const std::vector<char> *mask)
{
    comm = cm;
    B = b;
    BB = B*B;
    nx0 = nx;
    ny0 = ny;

    for(int q=0; q<4; ++q)
    nb[q] = n[q]>=0 ? n[q] : MPI_PROC_NULL;

    alone = true;
    for(int q=0; q<4; ++q)
    if(nb[q]!=MPI_PROC_NULL)
    alone = false;

    // the neighbours across a rank face must have the same number of columns along it (Cartesian partition)
    {
    int mine[2] = {ny, nx}, other[4] = {-1,-1,-1,-1};
    MPI_Sendrecv(&mine[0],1,MPI_INT,nb[0],7600,&other[1],1,MPI_INT,nb[1],7600,comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&mine[0],1,MPI_INT,nb[1],7601,&other[0],1,MPI_INT,nb[0],7601,comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&mine[1],1,MPI_INT,nb[2],7602,&other[3],1,MPI_INT,nb[3],7602,comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&mine[1],1,MPI_INT,nb[3],7603,&other[2],1,MPI_INT,nb[2],7603,comm,MPI_STATUS_IGNORE);

    int bad = 0;
    for(int q=0; q<2; ++q)
    if(nb[q]!=MPI_PROC_NULL && other[q]!=ny) bad = 1;
    for(int q=2; q<4; ++q)
    if(nb[q]!=MPI_PROC_NULL && other[q]!=nx) bad = 1;

    int gbad = 0;
    MPI_Allreduce(&bad,&gbad,1,MPI_INT,MPI_MAX,comm);

    if(gbad)
    {
    int rank; MPI_Comm_rank(comm,&rank);
    if(rank==0)
    std::cout<<"!!! NHFLOW A 510 11/12/13 with A 520 > 0: the multigrid of the coupled system needs a Cartesian partition (the same columns across rank faces) !!!"<<std::endl;
    exit(1);
    }
    }

    layout(mask);
    A0.assign((size_t)nx*ny*CS,0.0);
    xh0.assign((size_t)(nx+2)*(ny+2)*B,0.0);
    nz0.assign((size_t)5*BB,0);
    cf.reset();
    cd.reset();

    // the number of levels (the cores repeat it with their setup at the first finalize)
    int nl = 1, ax = nx, ay = ny;
    while((ax>2 || ay>2) && nl<30)
    {
    ax = (ax+1)/2;
    ay = (ay+1)/2;
    ++nl;
    }

    MPI_Allreduce(&nl,&nlev,1,MPI_INT,MPI_MIN,comm);
}

int nhflow_si_nhmg::levels() const
{
    return nlev;
}

void nhflow_si_nhmg::zero()
{
    std::fill(A0.begin(),A0.end(),0.0);
    std::fill(nz0.begin(),nz0.end(),0);
}

void nhflow_si_nhmg::mark_all()
{
    make_dense();
    std::fill(nz0.begin(),nz0.end(),1);
}

void nhflow_si_nhmg::layout(const std::vector<char> *mask)
{
    // the diagonal block in full, the others by the mask (all without one); the entries of a block in row-major order
    dense = (mask==nullptr);
    slot.assign((size_t)5*BB,-1);
    sl.assign(5,std::vector<std::pair<int,int>>());
    CS = 0;

    for(int d=0; d<5; ++d)
    {
        offs[d] = CS;
        int n = 0;

        for(int e=0; e<BB; ++e)
        if(d==0 || dense || (*mask)[(size_t)d*BB + e])
        {
        slot[(size_t)d*BB + e] = n;
        sl[d].push_back({e,n});
        ++n;
        }

        CS += n;
    }

    // the entries of the blocks outside the stored pattern (set_column)
    offp.assign(5,std::vector<int>());
    for(int d=0; d<5; ++d)
    for(int e=0; e<BB; ++e)
    if(slot[(size_t)d*BB + e]<0)
    offp[d].push_back(e);

    // the arrow layout (bmul): the neighbour blocks hold exactly row 0 and c = 0, r, r+1 of the rows r >= 1, in this order
    // (the mask of nhflow_momentum_SI::nh_mask). The products in the same order as with the pattern lists: the same
    // results bit for bit, 2D standing wave 64^2 x 4: the multigrid solve 1.2 s -> 0.8 s. Environment NHMG_GENERIC: the
    // pattern lists (check)
    arrow = !dense && getenv("NHMG_GENERIC")==nullptr;
    for(int d=1; d<5 && arrow; ++d)
    {
        std::vector<std::pair<int,int>> ex;
        int n = 0;
        for(int r=0; r<B; ++r)
        for(int c=0; c<B; ++c)
        if(r==0 || c==0 || c==r || c==r+1)
        ex.push_back({r*B+c,n++});
        arrow = (ex==sl[d]);
    }
}

void nhflow_si_nhmg::make_dense()
{
    // an entry outside the mask (or block() wanted): every block in full, the values kept; the preconditioner new
    if(dense)
    return;

    const std::vector<double> old = A0;
    const std::vector<int> os = slot;
    const int ocs = CS;
    int ooffs[5];
    for(int d=0; d<5; ++d) ooffs[d] = offs[d];

    layout(nullptr);
    A0.assign((size_t)nx0*ny0*CS,0.0);

    for(size_t q=0; q<(size_t)nx0*ny0; ++q)
    for(int d=0; d<5; ++d)
    for(int e=0; e<BB; ++e)
    if(os[(size_t)d*BB + e]>=0)
    A0[q*CS + offs[d] + slot[(size_t)d*BB + e]] = old[q*ocs + ooffs[d] + os[(size_t)d*BB + e]];

    cf.reset();
    cd.reset();
}

void nhflow_si_nhmg::set_column(int ix, int iy, int nd, const double *A, int r0, double s)
{
    // an entry outside the stored pattern that is not zero: every block in full (make_dense), then as set()
    if(!dense)
    for(int d=0; d<nd; ++d)
    for(const int e : offp[d])
    if(A[(size_t)d*BB + e]!=0.0)
    {
    make_dense();
    d = nd;
    break;
    }

    double *col = &A0[(size_t)(ix*ny0 + iy)*CS];
    char *nz = nz0.data();

    for(int d=0; d<nd; ++d)
    {
        const double *Ad = A + (size_t)d*BB;
        double *cd = col + offs[d];
        char *nzd = nz + (size_t)d*BB;

        for(const auto &e : sl[d])
        {
            const double v = Ad[e.first];

            if(v!=0.0)
            {
            cd[e.second] = (e.first>=r0*B) ? v*s : v;
            nzd[e.first] = 1;
            }
        }
    }
}

void nhflow_si_nhmg::scale_rows(int ix, int iy, int d, int r0, double s)
{
    double *a = &A0[(size_t)(ix*ny0 + iy)*CS + offs[d]];

    for(int r=r0; r<B; ++r)
    for(int c=0; c<B; ++c)
    {
    const int sl0 = slot[(size_t)d*BB + r*B + c];
    if(sl0>=0)
    a[sl0] *= s;
    }
}

void nhflow_si_nhmg::pattern()
{
    pat0.assign((size_t)5*B,std::vector<nhmg_pent>());

    for(int d=0; d<5; ++d)
    for(int r=0; r<B; ++r)
    for(int c=0; c<B; ++c)
    if(nz0[(size_t)d*BB + r*B + c])
    pat0[(size_t)d*B + r].push_back({c,slot[(size_t)d*BB + r*B + c]});
}

void nhflow_si_nhmg::finalize(bool refresh)
{
    pattern();

    if(single && !cf)
    {
    cf = std::make_shared<nhflow_si_nhmg_core<float>>();
    cf->setup(comm,B,nx0,ny0,nb,alone,this);
    refresh = true;
    }

    if(!single && !cd)
    {
    cd = std::make_shared<nhflow_si_nhmg_core<double>>();
    cd->setup(comm,B,nx0,ny0,nb,alone,this);
    refresh = true;
    }

    if(single)
    {
    cd.reset();      // the double levels only when used
    if(refresh) cf->build(A0,pat0);
    }

    else
    {
    cf.reset();
    if(refresh) cd->build(A0,pat0);
    }
}

void nhflow_si_nhmg::matvec(const double *x, double *y)
{
    for(int ix=0; ix<nx0; ++ix)
    for(int iy=0; iy<ny0; ++iy)
    memcpy(&xh0[(size_t)hix(ny0,ix,iy)*B],&x[(size_t)(ix*ny0 + iy)*B],B*sizeof(double));

    halo<double>(comm,nb,alone,B,nx0,ny0,xh0,sbuf,rbuf);

    dispatchB(B,[&]<int NB>()
    {
    for(int ix=0; ix<nx0; ++ix)
    for(int iy=0; iy<ny0; ++iy)
    {
        double *yy = &y[(size_t)(ix*ny0 + iy)*B];
        std::fill(yy,yy+B,0.0);

        for(int d=0; d<5; ++d)
        {
            const double *a = &A0[(size_t)(ix*ny0 + iy)*CS + offs[d]];
            const double *xx = &xh0[(size_t)hix(ny0,ix+dix[d],iy+diy[d])*B];
            bmulB<NB,double,double>(B,arrow,d,&pat0[(size_t)d*B],a,xx,yy,1.0);
        }
    }
    });
}

void nhflow_si_nhmg::vcycle(const double *r, double *z)
{
    if(single)
    {
    cf->o = this;
    cf->vcycle(r,z);
    }

    else
    {
    cd->o = this;
    cd->vcycle(r,z);
    }
}
