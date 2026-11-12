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

// y[r] += sign sum_c a[r B + c] x[c] over the nonzero pattern pat[r] of a block
template<class T, class U> inline void block_mul(int B, const std::vector<int> *pat, const T *a, const U *x, U *y, U sign)
{
    for(int r=0; r<B; ++r)
    {
        U t = 0;
        for(int c : pat[r])
        t += (U)a[r*B+c]*x[c];
        y[r] += sign*t;
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
        std::vector<T> A;       // 5 blocks per column
        std::vector<T> Di;      // the inverse of the diagonal block
        std::vector<T> x, b, r;     // with the halo of one column
        std::vector<T> kb, v1, w1, v2, w2;     // K-cycle
        std::vector<std::vector<int>> pat;      // the columns c of row r of block d that can be nonzero: pat[d*B + r]
    };

    MPI_Comm comm;
    int B=0, BB=0;
    int nb[4];
    bool alone = true;
    std::vector<level> L;
    std::vector<T> sbuf, rbuf, tmp, tmp2;
    const nhflow_si_nhmg *o = nullptr;     // the options

    void setup(MPI_Comm cm, int b, int nx, int ny, const int n[4], bool aln)
    {
        comm = cm; B = b; BB = B*B; alone = aln;
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
        v.A.assign((size_t)ax*ay*5*BB,T(0));
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
    void build(const std::vector<double> &A0, const std::vector<std::vector<int>> &pat0)
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

            for(size_t q=0; q<nc*5*BB; ++q)
            f.A[q] = (T)fine[q];

            for(size_t q=0; q<nc; ++q)
            {
            invert(B,&fine[q*5*BB],inv.data(),sc.data(),pv.data());
            for(int e=0; e<BB; ++e)
            f.Di[q*BB+e] = (T)inv[e];
            }

            if(l+1==(int)L.size())
            break;

            // the coarse operator: the blocks of the fine columns of an aggregate summed, a coupling within the aggregate
            // into its diagonal block; its pattern: block d (> 0) of the fine one, the diagonal block the union
            level &c = L[l+1];
            std::vector<double> &coarse = cdv;
            coarse.assign((size_t)c.nx*c.ny*5*BB,0.0);

            for(int ix=0; ix<f.nx; ++ix)
            for(int iy=0; iy<f.ny; ++iy)
            {
                const int Ix = ix/2, Iy = iy/2;

                for(int d=0; d<5; ++d)
                {
                    const double *a = &fine[((size_t)(ix*f.ny + iy)*5 + d)*BB];
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

                    double *t = &coarse[((size_t)(Ix*c.ny + Iy)*5 + dc)*BB];
                    for(int q=0; q<BB; ++q)
                    t[q] += a[q];
                }
            }

            c.pat = f.pat;
            for(int r=0; r<B; ++r)
            {
                std::vector<char> u(B,0);
                for(int d=0; d<5; ++d)
                for(int cc : f.pat[(size_t)d*B + r])
                u[cc] = 1;

                std::vector<int> &p0 = c.pat[r];
                p0.clear();
                for(int cc=0; cc<B; ++cc)
                if(u[cc]) p0.push_back(cc);
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

        for(int s=0; s<sweeps; ++s)
        for(int colour=0; colour<2; ++colour)
        {
            const bool first = zero && s==0 && colour==0;

            if(!first)
            exchange(v,v.x);

            for(int ix=0; ix<v.nx; ++ix)
            for(int iy=0; iy<v.ny; ++iy)
            {
                if(((ix+iy)&1)!=colour)
                continue;

                const size_t q = (size_t)(ix*v.ny + iy);
                const T *bb = &v.b[(size_t)hix(v.ny,ix,iy)*B];

                for(int r=0; r<B; ++r)
                tmp[r] = bb[r];

                if(!first)
                for(int d=1; d<5; ++d)
                {
                    const T *a = &v.A[(q*5 + d)*BB];
                    const T *xx = &v.x[(size_t)hix(v.ny,ix+dix[d],iy+diy[d])*B];
                    block_mul<T,T>(B,&v.pat[(size_t)d*B],a,xx,tmp.data(),T(-1));
                }

                const T *di = &v.Di[q*BB];
                T *xx = &v.x[(size_t)hix(v.ny,ix,iy)*B];
                for(int r=0; r<B; ++r)
                {
                    T t = 0;
                    for(int c=0; c<B; ++c)
                    t += di[r*B+c]*tmp[c];
                    xx[r] = t;
                }
            }
        }
    }

    // r = b - A x; skip >= 0: the colour relaxed last (its columns solved with the current neighbours: r = 0 there), the
    // columns at a rank face computed in full (the colours of the ranks need not match there)
    void residual(int lev, int skip=-1)
    {
        level &v = L[lev];
        exchange(v,v.x);

        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            const size_t q = (size_t)(ix*v.ny + iy);
            T *rr = &v.r[(size_t)hix(v.ny,ix,iy)*B];

            if(skip>=0 && ((ix+iy)&1)==skip && (alone || (ix>0 && ix<v.nx-1 && iy>0 && iy<v.ny-1)))
            {
            std::fill(rr,rr+B,T(0));
            continue;
            }

            const T *bb = &v.b[(size_t)hix(v.ny,ix,iy)*B];

            for(int r=0; r<B; ++r)
            rr[r] = bb[r];

            for(int d=0; d<5; ++d)
            {
                const T *a = &v.A[(q*5 + d)*BB];
                const T *xx = &v.x[(size_t)hix(v.ny,ix+dix[d],iy+diy[d])*B];
                block_mul<T,T>(B,&v.pat[(size_t)d*B],a,xx,rr,T(-1));
            }
        }
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
        exchange(v,x);

        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            const size_t q = (size_t)(ix*v.ny + iy);
            T *yy = &y[(size_t)hix(v.ny,ix,iy)*B];
            std::fill(yy,yy+B,T(0));

            for(int d=0; d<5; ++d)
            {
                const T *a = &v.A[(q*5 + d)*BB];
                const T *xx = &x[(size_t)hix(v.ny,ix+dix[d],iy+diy[d])*B];
                block_mul<T,T>(B,&v.pat[(size_t)d*B],a,xx,yy,T(1));
            }
        }
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
void nhflow_si_nhmg::setup(MPI_Comm cm, int b, int nx, int ny, const int n[4])
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

    A0.assign((size_t)nx*ny*5*BB,0.0);
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
    std::fill(nz0.begin(),nz0.end(),1);
}

void nhflow_si_nhmg::pattern()
{
    pat0.assign((size_t)5*B,std::vector<int>());

    for(int d=0; d<5; ++d)
    for(int r=0; r<B; ++r)
    for(int c=0; c<B; ++c)
    if(nz0[(size_t)d*BB + r*B + c])
    pat0[(size_t)d*B + r].push_back(c);
}

void nhflow_si_nhmg::finalize(bool refresh)
{
    pattern();

    if(single && !cf)
    {
    cf = std::make_shared<nhflow_si_nhmg_core<float>>();
    cf->setup(comm,B,nx0,ny0,nb,alone);
    refresh = true;
    }

    if(!single && !cd)
    {
    cd = std::make_shared<nhflow_si_nhmg_core<double>>();
    cd->setup(comm,B,nx0,ny0,nb,alone);
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

    for(int ix=0; ix<nx0; ++ix)
    for(int iy=0; iy<ny0; ++iy)
    {
        double *yy = &y[(size_t)(ix*ny0 + iy)*B];
        std::fill(yy,yy+B,0.0);

        for(int d=0; d<5; ++d)
        {
            const double *a = &A0[((size_t)(ix*ny0 + iy)*5 + d)*BB];
            const double *xx = &xh0[(size_t)hix(ny0,ix+dix[d],iy+diy[d])*B];
            block_mul<double,double>(B,&pat0[(size_t)d*B],a,xx,yy,1.0);
        }
    }
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
