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
}

void nhflow_si_nhmg::setup(MPI_Comm cm, int b, int nx, int ny, const int n[4])
{
    comm = cm;
    B = b;
    BB = B*B;

    for(int q=0; q<4; ++q)
    nb[q] = n[q]>=0 ? n[q] : MPI_PROC_NULL;

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

    // levels: 2 x 2 aggregates while a rank holds more than 2 columns in a direction; the same number on all ranks
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
    v.A.assign((size_t)ax*ay*5*BB,0.0);
    v.LU.assign((size_t)ax*ay*BB,0.0);
    v.piv.assign((size_t)ax*ay*B,0);
    v.x.assign((size_t)(ax+2)*(ay+2)*B,0.0);
    v.b = v.x;
    v.r = v.x;
    v.kb = v.v1 = v.w1 = v.v2 = v.w2 = v.x;
    ax = (ax+1)/2;
    ay = (ay+1)/2;
    }

    tmp.assign(B,0.0);
}

void nhflow_si_nhmg::zero()
{
    std::fill(L[0].A.begin(),L[0].A.end(),0.0);
}

void nhflow_si_nhmg::lu(double *a, int *pv)
{
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
}

void nhflow_si_nhmg::lusolve(const double *a, const int *pv, double *x)
{
    // the row swaps of the factorisation apply to the whole rows (L included): all of them first, then L y = P x
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
}

void nhflow_si_nhmg::finalize()
{
    // the Galerkin operators of the coarse levels: the blocks of the fine columns of an aggregate summed, a coupling
    // within the aggregate into its diagonal block
    for(int l=0; l+1<(int)L.size(); ++l)
    {
        level &f = L[l], &c = L[l+1];
        std::fill(c.A.begin(),c.A.end(),0.0);

        for(int ix=0; ix<f.nx; ++ix)
        for(int iy=0; iy<f.ny; ++iy)
        {
            const int Ix = ix/2, Iy = iy/2;

            for(int d=0; d<5; ++d)
            {
                const double *a = &f.A[((size_t)(ix*f.ny + iy)*5 + d)*BB];
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

                double *t = &c.A[((size_t)(Ix*c.ny + Iy)*5 + dc)*BB];
                for(int q=0; q<BB; ++q)
                t[q] += a[q];
            }
        }
    }

    // the LU of the diagonal blocks
    for(level &v : L)
    for(int q=0; q<v.nx*v.ny; ++q)
    {
    double *a = &v.LU[(size_t)q*BB];
    memcpy(a,&v.A[(size_t)q*5*BB],BB*sizeof(double));
    lu(a,&v.piv[(size_t)q*B]);
    }
}

void nhflow_si_nhmg::exchange(level &l, std::vector<double> &v)
{
    // x: the strips ix = 0 and nx-1 (ny columns), y: iy = 0 and ny-1 (nx columns)
    const int nx = l.nx, ny = l.ny;

    auto strip = [&](bool xdir, int idx, bool get, std::vector<double> &buf)
    {
        const int n = xdir ? ny : nx;
        buf.resize((size_t)n*B);

        for(int m=0; m<n; ++m)
        {
            const int ix = xdir ? idx : m, iy = xdir ? m : idx;
            double *p = &v[(size_t)hix(l,ix,iy)*B];

            if(get)
            memcpy(&buf[(size_t)m*B],p,B*sizeof(double));
            else
            memcpy(p,&buf[(size_t)m*B],B*sizeof(double));
        }
    };

    auto pass = [&](bool xdir, int send_idx, int dest, int recv_idx, int src, int tag)
    {
        const int n = (xdir ? ny : nx)*B;
        strip(xdir,send_idx,true,sbuf);
        rbuf.assign(n,0.0);
        MPI_Sendrecv(sbuf.data(),n,MPI_DOUBLE,dest,tag,rbuf.data(),n,MPI_DOUBLE,src,tag,comm,MPI_STATUS_IGNORE);

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

void nhflow_si_nhmg::matvec(const double *x, double *y)
{
    level &v = L[0];

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    memcpy(&v.x[(size_t)hix(v,ix,iy)*B],&x[(size_t)(ix*v.ny + iy)*B],B*sizeof(double));

    exchange(v,v.x);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    {
        double *yy = &y[(size_t)(ix*v.ny + iy)*B];
        std::fill(yy,yy+B,0.0);

        for(int d=0; d<5; ++d)
        {
            const double *a = &v.A[((size_t)(ix*v.ny + iy)*5 + d)*BB];
            const double *xx = &v.x[(size_t)hix(v,ix+dix[d],iy+diy[d])*B];

            for(int r=0; r<B; ++r)
            {
            double s = 0.0;
            for(int c=0; c<B; ++c)
            s += a[r*B+c]*xx[c];
            yy[r] += s;
            }
        }
    }
}

void nhflow_si_nhmg::smooth(int lev, int sweeps)
{
    level &v = L[lev];

    for(int s=0; s<sweeps; ++s)
    for(int colour=0; colour<2; ++colour)
    {
        exchange(v,v.x);

        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            if(((ix+iy)&1)!=colour)
            continue;

            const size_t q = (size_t)(ix*v.ny + iy);
            const double *bb = &v.b[(size_t)hix(v,ix,iy)*B];

            for(int r=0; r<B; ++r)
            tmp[r] = bb[r];

            for(int d=1; d<5; ++d)
            {
                const double *a = &v.A[(q*5 + d)*BB];
                const double *xx = &v.x[(size_t)hix(v,ix+dix[d],iy+diy[d])*B];

                for(int r=0; r<B; ++r)
                {
                double t = 0.0;
                for(int c=0; c<B; ++c)
                t += a[r*B+c]*xx[c];
                tmp[r] -= t;
                }
            }

            lusolve(&v.LU[q*BB],&v.piv[q*B],tmp.data());
            memcpy(&v.x[(size_t)hix(v,ix,iy)*B],tmp.data(),B*sizeof(double));
        }
    }
}

void nhflow_si_nhmg::residual(int lev)
{
    level &v = L[lev];
    exchange(v,v.x);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    {
        const size_t q = (size_t)(ix*v.ny + iy);
        double *rr = &v.r[(size_t)hix(v,ix,iy)*B];
        const double *bb = &v.b[(size_t)hix(v,ix,iy)*B];

        for(int r=0; r<B; ++r)
        rr[r] = bb[r];

        for(int d=0; d<5; ++d)
        {
            const double *a = &v.A[(q*5 + d)*BB];
            const double *xx = &v.x[(size_t)hix(v,ix+dix[d],iy+diy[d])*B];

            for(int r=0; r<B; ++r)
            {
            double t = 0.0;
            for(int c=0; c<B; ++c)
            t += a[r*B+c]*xx[c];
            rr[r] -= t;
            }
        }
    }
}

void nhflow_si_nhmg::cycle(int lev)
{
    level &v = L[lev];
    std::fill(v.x.begin(),v.x.end(),0.0);

    if(lev==(int)L.size()-1)
    {
    smooth(lev,ncoarse);
    return;
    }

    smooth(lev,npre);
    residual(lev);

    level &c = L[lev+1];
    std::fill(c.b.begin(),c.b.end(),0.0);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    {
        double *cb = &c.b[(size_t)hix(c,ix/2,iy/2)*B];
        const double *rr = &v.r[(size_t)hix(v,ix,iy)*B];
        for(int r=0; r<B; ++r)
        cb[r] += rr[r];
    }

    if(kcycle && lev+1<(int)L.size()-1)
    kcorrect(lev+1);
    else
    cycle(lev+1);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    {
        double *xx = &v.x[(size_t)hix(v,ix,iy)*B];
        const double *cx = &c.x[(size_t)hix(c,ix/2,iy/2)*B];
        for(int r=0; r<B; ++r)
        xx[r] += omega_c*cx[r];
    }

    smooth(lev,npost);
}

void nhflow_si_nhmg::apply(int lev, std::vector<double> &x, std::vector<double> &y)
{
    level &v = L[lev];
    exchange(v,x);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    {
        const size_t q = (size_t)(ix*v.ny + iy);
        double *yy = &y[(size_t)hix(v,ix,iy)*B];
        std::fill(yy,yy+B,0.0);

        for(int d=0; d<5; ++d)
        {
            const double *a = &v.A[(q*5 + d)*BB];
            const double *xx = &x[(size_t)hix(v,ix+dix[d],iy+diy[d])*B];

            for(int r=0; r<B; ++r)
            {
            double t = 0.0;
            for(int c=0; c<B; ++c)
            t += a[r*B+c]*xx[c];
            yy[r] += t;
            }
        }
    }
}

void nhflow_si_nhmg::dots(int lev, const std::vector<double> *const *u, const std::vector<double> *const *w, int n, double *res)
{
    level &v = L[lev];

    for(int k=0; k<n; ++k)
    {
        double s = 0.0;
        for(int ix=0; ix<v.nx; ++ix)
        for(int iy=0; iy<v.ny; ++iy)
        {
            const size_t h = (size_t)hix(v,ix,iy)*B;
            for(int r=0; r<B; ++r)
            s += (*u[k])[h+r]*(*w[k])[h+r];
        }
        res[k] = s;
    }

    MPI_Allreduce(MPI_IN_PLACE,res,n,MPI_DOUBLE,MPI_SUM,comm);
}

// K-cycle (Notay): the coarse-level correction by two steps of GCR on A_lev e = b_lev, preconditioned by the cycle of
// the level; the second step only when the first leaves more than a quarter of the residual
void nhflow_si_nhmg::kcorrect(int lev)
{
    level &v = L[lev];
    v.kb = v.b;

    cycle(lev);
    v.v1 = v.x;
    apply(lev,v.v1,v.w1);

    double d[3];
    {
    const std::vector<double> *u[3] = {&v.w1,&v.w1,&v.kb};
    const std::vector<double> *w[3] = {&v.w1,&v.kb,&v.kb};
    dots(lev,u,w,3,d);
    }

    const double rho1 = d[0];
    if(rho1<=0.0)
    return;

    const double a1 = d[1]/rho1;
    const double rn2 = d[2] - d[1]*d[1]/rho1;     // |b - a1 w1|^2

    if(rn2 <= 0.0625*d[2])
    {
        for(size_t q=0; q<v.x.size(); ++q)
        v.x[q] = a1*v.v1[q];
        return;
    }

    for(size_t q=0; q<v.b.size(); ++q)
    v.b[q] = v.kb[q] - a1*v.w1[q];

    cycle(lev);
    v.v2 = v.x;
    apply(lev,v.v2,v.w2);

    {
    const std::vector<double> *u[3] = {&v.w1,&v.w2,&v.w2};
    const std::vector<double> *w[3] = {&v.w2,&v.w2,&v.b};
    dots(lev,u,w,3,d);
    }

    const double g = d[0]/rho1;
    const double rho2 = d[1] - d[0]*d[0]/rho1;            // |w2 - g w1|^2
    const double a2 = rho2>1.0e-300*d[1] ? (d[2] - g*0.0)/rho2 : 0.0;    // (w2 - g w1).r2 = w2.r2 (r2 orthogonal to w1)

    for(size_t q=0; q<v.x.size(); ++q)
    v.x[q] = a1*v.v1[q] + a2*(v.v2[q] - g*v.v1[q]);
}

void nhflow_si_nhmg::vcycle(const double *r, double *z)
{
    level &v = L[0];
    std::fill(v.b.begin(),v.b.end(),0.0);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    memcpy(&v.b[(size_t)hix(v,ix,iy)*B],&r[(size_t)(ix*v.ny + iy)*B],B*sizeof(double));

    cycle(0);

    for(int ix=0; ix<v.nx; ++ix)
    for(int iy=0; iy<v.ny; ++iy)
    memcpy(&z[(size_t)(ix*v.ny + iy)*B],&v.x[(size_t)hix(v,ix,iy)*B],B*sizeof(double));
}
