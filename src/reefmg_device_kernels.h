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

#ifndef REEFMG_DEVICE_KERNELS_H_
#define REEFMG_DEVICE_KERNELS_H_

//  The kernels of the REEFMG device path (N 16), written once and compiled
//  twice:
//
//    - as C++ in namespace rmg_kern, for the reference device (N 16 2) that
//      runs them in plain loops on the CPU, and
//    - as Metal Shading Language (MSL) for the Apple GPU (N 16 1): the same
//      text, stringified into rmg_kernel_source and compiled at run time by
//      the Metal framework (reefmg_device_metal.mm adds the entry points).
//      No offline Metal compiler is needed.
//
//  The text must stay inside the common subset of C++ and MSL: int,
//  unsigned, float, structs of int, no long, no double (Apple GPUs have no
//  fp64), no library calls, no preprocessor directives, no string literals.
//  Pointers to level arrays carry the MSL address space 'device'; in C++ the
//  word is defined away while the text is expanded.  The stringified copy is
//  taken by the # operator, before expansion, so the MSL keeps it.
//
//  Device layout (not reefmg_core's): k SLOWEST, one halo cell in x and y,
//  none in z,
//      q = k*sk + (i+1)*si + (j+1),   si = ny+2,  sk = (nx+2)*(ny+2),
//  so that neighbouring work items (neighbouring j) touch neighbouring
//  memory.  A level holds fewer than 2^31 cells.
//
//  One work item each:
//    rmg_smooth    one column of one colour ((i+j)&1): red-black (zebra) line
//                  Gauss-Seidel with the stored column factorisation, in
//                  place; columns without an active cell are skipped
//    rmg_residual  one interior cell: r = f - A u
//    rmg_restrict  one interior coarse cell: f = mean of the fine residuals
//                  of its children, u = 0
//    rmg_prolong   one interior fine cell: u += P u_coarse, masked by act
//    rmg_coarsen   one interior coarse cell: the coarse operator from the
//                  fine one (reefmg_core::coarsen_t) and its activity
//    rmg_factor    one column: the LU factors of its tridiagonal matrix
//                  (reefmg_core::factor_lines_t)
//    rmg_wrap      one halo cell of a periodic direction
//    rmg_zero      one entry of an array

#define RMG_DEVICE_SOURCE(...)                                               \
    static const char *const rmg_kernel_source = #__VA_ARGS__;               \
    namespace rmg_kern { __VA_ARGS__ }

#define device

RMG_DEVICE_SOURCE(

//  one level
struct rmg_lp
{
    int nx, ny, nz;
    int si, sk;
    int perx, pery;
};

//  a fine level and the next coarser one
struct rmg_tp
{
    int fnx, fny, nz, fsi, fsk;
    int cnx, cny, csi, csk;
    int rx, ry;
    int perx, pery;
};

inline int rmg_q(int i, int j, int k, int si, int sk)
{
    return k*sk + (i+1)*si + (j+1);
}

//  columns per colour: work items of rmg_smooth
inline int rmg_ncolour(const rmg_lp P)
{
    return P.nx*((P.ny+1)/2);
}

inline void rmg_smooth(const rmg_lp P, int colour, unsigned t,
                       const device float *f, device float *u,
                       const device float *cn, const device float *cs,
                       const device float *cw, const device float *ce,
                       const device float *cb, const device float *ctc,
                       const device float *cti, const device unsigned char *act)
{
    const int nyh = (P.ny+1)/2;
    const int it = int(t);
    const int i = it / nyh;
    if(i>=P.nx)
    return;
    const int j = ((i+colour)&1) + 2*(it % nyh);
    if(j>=P.ny)
    return;

    const int si = P.si, sk = P.sk, nz = P.nz;
    const int q0 = (i+1)*si + (j+1);

    bool any = false;
    for(int k=0; k<nz; ++k)
    if(act[q0+k*sk]!=0)
    {
        any = true;
        break;
    }
    if(!any)
    return;

    float prev = 0.0f;
    for(int k=0; k<nz; ++k)
    {
        const int q = q0 + k*sk;
        float r = f[q] - cn[q]*u[q+si] - cs[q]*u[q-si] - cw[q]*u[q+1] - ce[q]*u[q-1];
        if(k>0)
        r -= cb[q]*prev;
        r *= cti[q];
        u[q] = r;
        prev = r;
    }

    for(int k=nz-2; k>=0; --k)
    {
        const int q = q0 + k*sk;
        u[q] -= ctc[q]*u[q+sk];
    }
}

inline void rmg_residual(const rmg_lp P, unsigned t,
                         const device float *f, const device float *u, device float *r,
                         const device float *cp, const device float *cn,
                         const device float *cs, const device float *cw,
                         const device float *ce, const device float *ct,
                         const device float *cb)
{
    const int it = int(t);
    const int j = it % P.ny;
    const int c = it / P.ny;
    const int i = c % P.nx;
    const int k = c / P.nx;
    const int si = P.si, sk = P.sk;
    const int q = rmg_q(i,j,k,si,sk);

    float y = cp[q]*u[q] + cn[q]*u[q+si] + cs[q]*u[q-si] + cw[q]*u[q+1] + ce[q]*u[q-1];
    if(k<P.nz-1)
    y += ct[q]*u[q+sk];
    if(k>0)
    y += cb[q]*u[q-sk];

    r[q] = f[q] - y;
}

//  children summed in the order of reefmg_core::restrict_xy:
//  (a,b) = (0,0),(0,1),(1,0),(1,1)
inline void rmg_restrict(const rmg_tp T, unsigned t,
                         const device float *fr, device float *cf, device float *cu)
{
    const int it = int(t);
    const int J = it % T.cny;
    const int c = it / T.cny;
    const int I = c % T.cnx;
    const int k = c / T.cnx;

    const int i0 = T.rx*I, j0 = T.ry*J;
    const int i1 = (T.rx==2 && i0+1<T.fnx) ? i0+1 : -1;
    const int j1 = (T.ry==2 && j0+1<T.fny) ? j0+1 : -1;

    const float r00 = fr[rmg_q(i0,j0,k,T.fsi,T.fsk)];
    float s;

    if(i1>=0 && j1>=0)
    s = (((r00 + fr[rmg_q(i0,j1,k,T.fsi,T.fsk)]) + fr[rmg_q(i1,j0,k,T.fsi,T.fsk)])
              + fr[rmg_q(i1,j1,k,T.fsi,T.fsk)])*0.25f;
    else if(i1>=0)
    s = (r00 + fr[rmg_q(i1,j0,k,T.fsi,T.fsk)])*0.5f;
    else if(j1>=0)
    s = (r00 + fr[rmg_q(i0,j1,k,T.fsi,T.fsk)])*0.5f;
    else
    s = r00;

    const int qc = rmg_q(I,J,k,T.csi,T.csk);
    cf[qc] = s;
    cu[qc] = 0.0f;
}

//  cell-centred bilinear in the plane, identity in z (reefmg_core::prolong_xy)
inline void rmg_prolong(const rmg_tp T, unsigned t,
                        const device float *cu, device float *fu,
                        const device unsigned char *act)
{
    const int it = int(t);
    const int j = it % T.fny;
    const int c = it / T.fny;
    const int i = c % T.fnx;
    const int k = c / T.fnx;

    const int q = rmg_q(i,j,k,T.fsi,T.fsk);
    if(act[q]==0)
    return;

    const int I = i / T.rx, J = j / T.ry;
    int Ia = I, Ja = J;
    float wx = 1.0f, wy = 1.0f;

    if(T.rx==2)
    {
        Ia = (i%2==0) ? I-1 : I+1;
        wx = 0.75f;
    }
    if(T.ry==2)
    {
        Ja = (j%2==0) ? J-1 : J+1;
        wy = 0.75f;
    }
    Ia = Ia<0 ? 0 : (Ia>T.cnx-1 ? T.cnx-1 : Ia);
    Ja = Ja<0 ? 0 : (Ja>T.cny-1 ? T.cny-1 : Ja);

    const float w00 = wx*wy, w10 = (1.0f-wx)*wy, w01 = wx*(1.0f-wy), w11 = (1.0f-wx)*(1.0f-wy);

    const float c00 = cu[rmg_q(I ,J ,k,T.csi,T.csk)];
    const float c10 = cu[rmg_q(Ia,J ,k,T.csi,T.csk)];
    const float c01 = cu[rmg_q(I ,Ja,k,T.csi,T.csk)];
    const float c11 = cu[rmg_q(Ia,Ja,k,T.csi,T.csk)];

    fu[q] += w00*c00 + w10*c10 + w01*c01 + w11*c11;
}

//  Coarse operator of one coarse cell, reefmg_core::coarsen_t in single
//  precision: the face coefficients recovered from the fine rows,
//  averaged over the fine faces of the coarse face and rebuilt with the
//  coarse widths; the row sum carried over.  Widths hx/hy have one ghost
//  entry each side, index i+1 for cell i.
inline void rmg_coarsen(const rmg_tp T, unsigned t,
                        const device float *fp, const device float *fn, const device float *fs,
                        const device float *fw, const device float *fe, const device float *ft,
                        const device float *fb, const device unsigned char *fact,
                        const device float *fhx, const device float *fhy,
                        device float *cp, device float *cn, device float *cs,
                        device float *cw, device float *ce, device float *ct,
                        device float *cb, device unsigned char *cact,
                        const device float *chx, const device float *chy)
{
    const int it = int(t);
    const int J = it % T.cny;
    const int c = it / T.cny;
    const int I = c % T.cnx;
    const int k = c / T.cnx;
    const int rx = T.rx, ry = T.ry;

    int ihi = rx*I, jhi = ry*J;
    for(int a=0; a<rx; ++a) if(rx*I+a<T.fnx) ihi = rx*I+a;
    for(int b=0; b<ry; ++b) if(ry*J+b<T.fny) jhi = ry*J+b;

    float volsum = 0.0f, tsum = 0.0f, bsum = 0.0f, xsum = 0.0f;
    float an = 0.0f, as = 0.0f, aw = 0.0f, ae = 0.0f;
    float wn = 0.0f, ws = 0.0f, ww = 0.0f, we = 0.0f;
    int na = 0;

    for(int i=rx*I; i<=ihi; ++i)
    for(int j=ry*J; j<=jhi; ++j)
    {
        const int qf = rmg_q(i,j,k,T.fsi,T.fsk);
        if(fact[qf]==0)
        continue;
        ++na;

        const float vp = fp[qf], vn = fn[qf], vs = fs[qf], vw = fw[qf];
        const float ve = fe[qf], vt = ft[qf], vb = fb[qf];
        const float hxf = fhx[i+1], hyf = fhy[j+1];
        const float vol = hxf*hyf;

        volsum += vol;
        tsum += vt*vol;
        bsum += vb*vol;
        xsum += (vp+vn+vs+vw+ve+vt+vb)*vol;

        if(i==ihi)
        {
            const float d = 0.5f*(hxf+fhx[i+2]);
            an += (-vn*d*hxf)*hyf;
            wn += hyf;
        }
        if(i==rx*I)
        {
            const float d = 0.5f*(hxf+fhx[i]);
            as += (-vs*d*hxf)*hyf;
            ws += hyf;
        }
        if(j==jhi)
        {
            const float d = 0.5f*(hyf+fhy[j+2]);
            aw += (-vw*d*hyf)*hxf;
            ww += hxf;
        }
        if(j==ry*J)
        {
            const float d = 0.5f*(hyf+fhy[j]);
            ae += (-ve*d*hyf)*hxf;
            we += hxf;
        }
    }

    const int qc = rmg_q(I,J,k,T.csi,T.csk);

    if(na==0)
    {
        cp[qc] = 1.0f;
        cn[qc] = 0.0f; cs[qc] = 0.0f; cw[qc] = 0.0f; ce[qc] = 0.0f; ct[qc] = 0.0f; cb[qc] = 0.0f;
        cact[qc] = 0;
        return;
    }
    cact[qc] = 1;

    const float HX = chx[I+1], HY = chy[J+1];
    const float DXn = 0.5f*(HX+chx[I+2]), DXs = 0.5f*(HX+chx[I]);
    const float DYw = 0.5f*(HY+chy[J+2]), DYe = 0.5f*(HY+chy[J]);

    float vt = tsum/volsum;
    float vb = bsum/volsum;
    float vn = (wn>0.0f) ? -(an/wn)/(DXn*HX) : 0.0f;
    float vs = (ws>0.0f) ? -(as/ws)/(DXs*HX) : 0.0f;
    float vw = (ww>0.0f) ? -(aw/ww)/(DYw*HY) : 0.0f;
    float ve = (we>0.0f) ? -(ae/we)/(DYe*HY) : 0.0f;

    if(I==0       && T.perx==0) vs = 0.0f;
    if(I==T.cnx-1 && T.perx==0) vn = 0.0f;
    if(J==0       && T.pery==0) ve = 0.0f;
    if(J==T.cny-1 && T.pery==0) vw = 0.0f;
    if(k==0)        vb = 0.0f;
    if(k==T.nz-1)   vt = 0.0f;

    cp[qc] = xsum/volsum - (vn+vs+vw+ve+vt+vb);
    cn[qc] = vn; cs[qc] = vs; cw[qc] = vw; ce[qc] = ve; ct[qc] = vt; cb[qc] = vb;
}

//  LU factors of one column, reefmg_core::factor_lines_t
inline void rmg_factor(const rmg_lp P, unsigned t,
                       const device float *cp, const device float *ct, const device float *cb,
                       device float *ctc, device float *cti, const device unsigned char *act)
{
    const int it = int(t);
    const int j = it % P.ny;
    const int i = it / P.ny;
    const int sk = P.sk;
    const int q0 = (i+1)*P.si + (j+1);

    bool any = false;
    for(int k=0; k<P.nz; ++k)
    if(act[q0+k*sk]!=0)
    {
        any = true;
        break;
    }
    if(!any)
    return;

    float inv = 1.0f/cp[q0];
    float tcp = ct[q0]*inv;
    cti[q0] = inv;
    ctc[q0] = tcp;

    for(int k=1; k<P.nz; ++k)
    {
        const int q = q0 + k*sk;
        inv = 1.0f/(cp[q] - cb[q]*tcp);
        tcp = ct[q]*inv;
        cti[q] = inv;
        ctc[q] = tcp;
    }
}

//  periodic halo: work items 0 .. ny*nz-1 the x direction, then nx*nz the y
//  direction; the halo corners are never read
inline void rmg_wrap(const rmg_lp P, unsigned t, device float *u)
{
    const int it = int(t);
    const int nx_ = P.nx, ny_ = P.ny, nz_ = P.nz;

    if(it < ny_*nz_)
    {
        if(P.perx==0)
        return;
        const int j = it % ny_;
        const int k = it / ny_;
        u[rmg_q(-1 ,j,k,P.si,P.sk)] = u[rmg_q(nx_-1,j,k,P.si,P.sk)];
        u[rmg_q(nx_,j,k,P.si,P.sk)] = u[rmg_q(0    ,j,k,P.si,P.sk)];
        return;
    }

    const int iy = it - ny_*nz_;
    if(P.pery==0 || iy>=nx_*nz_)
    return;
    const int i = iy % nx_;
    const int k = iy / nx_;
    u[rmg_q(i,-1 ,k,P.si,P.sk)] = u[rmg_q(i,ny_-1,k,P.si,P.sk)];
    u[rmg_q(i,ny_,k,P.si,P.sk)] = u[rmg_q(i,0    ,k,P.si,P.sk)];
}

inline int rmg_nwrap(const rmg_lp P)
{
    return P.ny*P.nz + P.nx*P.nz;
}

inline void rmg_zero(unsigned t, device float *a)
{
    a[t] = 0.0f;
}

)

#undef device

#endif
