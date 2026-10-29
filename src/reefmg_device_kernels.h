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
//      the Metal framework (reefmg_device_metal.cpp adds the entry points).
//      No offline Metal compiler is needed.
//
//  The text must stay inside the common subset of C++ and MSL: int,
//  unsigned, float, structs of int, no long, no double (Apple GPUs have no
//  fp64), no library calls, no preprocessor directives, no string literals.
//  Pointers to level arrays carry the MSL address space 'device'; in C++ the
//  word is defined away while the text is expanded (and 'constant', the
//  address space of the tail's parameters, becomes const).  The stringified
//  copy is taken by the # operator, before expansion, so the MSL keeps them.
//  Four names are defined outside the text, for C++ below and for MSL in
//  reefmg_device_metal.cpp: rmg_h (half precision storage: half in MSL, its
//  16 bits in C++), rmg_hf (rmg_h to float), rmg_toh (float to rmg_h, round
//  to nearest even) and rmg_barrier (a threadgroup barrier; nothing in C++,
//  where one work item does all the work of the group).
//
//  Device layout (not reefmg_core's): k SLOWEST, one halo cell in x and y,
//  none in z, and every row (i,k) split by the red-black colour (i+j)&1:
//
//      q = k*sk + (i+1)*si + c*hh + (j+1)/2,   c = (i+j)&1,
//      hh = (ny+3)/2 entries per colour,  si = 2*hh,  sk = (nx+2)*si.
//
//  A half-sweep of one colour then reads and writes one contiguous half of
//  every row - its own coefficients only - and its j-neighbours are the
//  contiguous other half: no cache line is loaded for the colour that is
//  not being smoothed.  The halo positions j = -1 and j = ny exist in both
//  layouts.  Work items run on a grid in memory order (the place in a half,
//  then colour and row, then the layer), see rmg_cell; a level holds fewer
//  than 2^31 entries.  rmg_q is also what reefmg_core uses to write the fine
//  level.
//
//  One work item each:
//    rmg_smooth    one column of one colour: red-black (zebra) line
//                  Gauss-Seidel with the stored column factorisation, in
//                  place; columns without an active cell are skipped; the
//                  first half-sweep of a cycle from u = 0 without reading u
//    rmg_smooth_h  the same with the coefficients in half precision (N 21,
//                  the default): every coupling divided by the pivot of the
//                  column's LU factorisation, dimensionless, so 16 bits hold
//                  it to 5e-4 without a range problem
//    rmg_restrict_residual
//                  one interior coarse cell: f = mean of the fine residuals
//                  f - A u of its children, computed in place, u = 0
//    rmg_prolong   one interior fine cell: u += P u_coarse, masked by act
//    rmg_coarsen   one interior coarse cell: the coarse operator from the
//                  fine one (reefmg_core::coarsen_t) and its activity
//    rmg_factor    one column: the LU factors of its tridiagonal matrix
//                  (reefmg_core::factor_lines_t)
//    rmg_factor_h  the same, and the scaled couplings of rmg_smooth_h
//    rmg_wrap      one halo cell of a periodic direction
//    rmg_zero      one entry of an array
//
//  And one threadgroup:
//    rmg_tail_run  the whole V-cycle of the small coarse levels at the bottom
//                  of the hierarchy - smoothing, residual-restriction,
//                  prolongation, the coarsest level - with barriers instead
//                  of one dispatch per step (reefmg_device::cycle)

//  C++ side of the four names
namespace rmg_kern
{
    typedef unsigned short rmg_h;

    inline float rmg_hf(rmg_h h)
    {
        union {unsigned u; float f;} o, magic;
        magic.u = 113u << 23;
        const unsigned shifted_exp = 0x7c00u << 13;
        o.u = (unsigned)(h & 0x7fffu) << 13;
        const unsigned e = shifted_exp & o.u;
        o.u += (127u - 15u) << 23;
        if(e == shifted_exp)            // inf, nan
        o.u += (128u - 16u) << 23;
        else if(e == 0)                 // zero, subnormal
        {
            o.u += 1u << 23;
            o.f -= magic.f;
        }
        o.u |= (unsigned)(h & 0x8000u) << 16;
        return o.f;
    }

    //  round to nearest even, as the GPU converts
    inline rmg_h rmg_toh(float x)
    {
        union {unsigned u; float f;} v, inf32, max16, dmagic;
        inf32.u = 255u << 23;
        max16.u = (127u + 16u) << 23;
        dmagic.u = ((127u - 15u) + (23u - 10u) + 1u) << 23;
        v.f = x;
        const unsigned sign = v.u & 0x80000000u;
        v.u ^= sign;
        unsigned o;
        if(v.u >= max16.u)
        o = (v.u > inf32.u) ? 0x7e00u : 0x7c00u;
        else if(v.u < (113u << 23))     // subnormal half
        {
            v.f += dmagic.f;
            o = v.u - dmagic.u;
        }
        else
        {
            const unsigned odd = (v.u >> 13) & 1u;
            v.u += ((unsigned)(15 - 127) << 23) + 0xfffu;
            v.u += odd;
            o = v.u >> 13;
        }
        return (rmg_h)(o | (sign >> 16));
    }

    inline void rmg_barrier(){}
}

#define RMG_DEVICE_SOURCE(...)                                               \
    static const char *const rmg_kernel_source = #__VA_ARGS__;               \
    namespace rmg_kern { __VA_ARGS__ }

#define device
#define constant const

RMG_DEVICE_SOURCE(

//  one level
struct rmg_lp
{
    int nx, ny, nz;
    int hh, si, sk;
    int perx, pery;
    int h16;          // the smoother's coefficients in half precision
};

//  a fine level and the next coarser one
struct rmg_tp
{
    int fnx, fny, nz, fhh, fsi, fsk;
    int cnx, cny, chh, csi, csk;
    int rx, ry;
    int perx, pery;
};

inline int rmg_q(int i, int j, int k, int hh, int si, int sk)
{
    return k*sk + (i+1)*si + ((i+j)&1)*hh + ((j+1)>>1);
}

//  Work items run on a grid (m, r, k): m the place in a colour half of a row
//  (0 .. hh-1), r = 2*i + c the row i and colour c (0 .. 2*nx-1), k the
//  layer - memory order, and no integer division.  rmg_cell gives the cell;
//  ok is false for the halo positions of a row, which have no work.
struct rmg_ijk
{
    int i, j, k, c, b;   // b = (j+1)&1: the parity of the place in the row
    bool ok;
};

inline rmg_ijk rmg_cell(int m, int r, int k, int ny)
{
    rmg_ijk x;
    x.i = r >> 1;
    x.c = r & 1;
    x.k = k;
    x.b = (x.c + x.i + 1) & 1;
    x.j = 2*m + x.b - 1;
    x.ok = (x.j>=0 && x.j<ny);
    return x;
}

//  r = f - A u at cell (i,j,k) of colour c: the j-neighbours are the other
//  colour half of the row, at places m+b (j+1) and m+b-1 (j-1); the
//  i-neighbours are the other half of the neighbouring rows, same place
inline float rmg_resid_at(const rmg_lp P, const rmg_ijk x,
                          const device float *f, const device float *u,
                          const device float *cp, const device float *cn,
                          const device float *cs, const device float *cw,
                          const device float *ce, const device float *ct,
                          const device float *cb)
{
    const int q = rmg_q(x.i,x.j,x.k,P.hh,P.si,P.sk);
    const int o = (1-2*x.c)*P.hh;
    float y = cp[q]*u[q] + cn[q]*u[q+P.si+o] + cs[q]*u[q-P.si+o]
                         + cw[q]*u[q+o+x.b] + ce[q]*u[q+o+x.b-1];
    if(x.k<P.nz-1)
    y += ct[q]*u[q+P.sk];
    if(x.k>0)
    y += cb[q]*u[q-P.sk];
    return f[q] - y;
}

//  Red-black line Gauss-Seidel, the column (i, place m) of one colour.
//  zero: the first half-sweep from u = 0 - the other colour is taken as zero
//  and set to zero (place m of the other half, every layer), which replaces
//  zeroing u before the cycle.
inline void rmg_smooth(const rmg_lp P, int colour, int m, int i, int zero,
                       const device float *f, device float *u,
                       const device float *cn, const device float *cs,
                       const device float *cw, const device float *ce,
                       const device float *cb, const device float *ctc,
                       const device float *cti, const device unsigned char *act)
{
    if(i>=P.nx || m>=P.hh)
    return;
    const int b = (colour+i+1)&1;
    const int j = 2*m + b - 1;

    const int si = P.si, sk = P.sk, nz = P.nz;
    const int q0 = (i+1)*si + colour*P.hh + m;
    const int o  = (1-2*colour)*P.hh;    // to the other colour of the row
    const int dw = o + b;                // j+1: place m+b of the other half
    const int de = o + b - 1;            // j-1: place m+b-1

    if(zero!=0)
    for(int k=0; k<nz; ++k)
    u[q0+o+k*sk] = 0.0f;

    if(j<0 || j>=P.ny)
    return;

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
        float r = f[q];
        if(zero==0)
        r -= cn[q]*u[q+si+o] + cs[q]*u[q-si+o] + cw[q]*u[q+dw] + ce[q]*u[q+de];
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

//  rmg_smooth with the factorisation folded into the couplings (rmg_factor_h):
//  an..ab = cn..cb times the inverse pivot ti of the cell, hc = ctc, all in
//  half precision; f is scaled by ti in single precision.  The same
//  forward and backward substitution as rmg_smooth, to the rounding of the
//  16-bit couplings (relative 2^-11).
inline void rmg_smooth_h(const rmg_lp P, int colour, int m, int i, int zero,
                         const device float *f, device float *u,
                         const device rmg_h *an, const device rmg_h *as,
                         const device rmg_h *aw, const device rmg_h *ae,
                         const device rmg_h *ab, const device rmg_h *hc,
                         const device float *cti, const device unsigned char *act)
{
    if(i>=P.nx || m>=P.hh)
    return;
    const int b = (colour+i+1)&1;
    const int j = 2*m + b - 1;

    const int si = P.si, sk = P.sk, nz = P.nz;
    const int q0 = (i+1)*si + colour*P.hh + m;
    const int o  = (1-2*colour)*P.hh;
    const int dw = o + b;
    const int de = o + b - 1;

    if(zero!=0)
    for(int k=0; k<nz; ++k)
    u[q0+o+k*sk] = 0.0f;

    if(j<0 || j>=P.ny)
    return;

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
        float r = f[q]*cti[q];
        if(zero==0)
        r -= rmg_hf(an[q])*u[q+si+o] + rmg_hf(as[q])*u[q-si+o] + rmg_hf(aw[q])*u[q+dw] + rmg_hf(ae[q])*u[q+de];
        if(k>0)
        r -= rmg_hf(ab[q])*prev;
        u[q] = r;
        prev = r;
    }

    for(int k=nz-2; k>=0; --k)
    {
        const int q = q0 + k*sk;
        u[q] -= rmg_hf(hc[q])*u[q+sk];
    }
}

//  Coarse cell (I,J,k): f = mean of the residuals f - A u of its fine
//  children, computed here (the fine residual is never stored), summed in
//  the order of reefmg_core::restrict_xy: (a,b) = (0,0),(0,1),(1,0),(1,1);
//  u = 0.
inline void rmg_restrict_residual(const rmg_tp T, const rmg_lp F, int M, int R, int K,
                                  const device float *f, const device float *u,
                                  const device float *cp, const device float *cn,
                                  const device float *cs, const device float *cw,
                                  const device float *ce, const device float *ct,
                                  const device float *cb,
                                  device float *cf, device float *cu)
{
    if(M>=T.chh)
    return;
    const rmg_ijk X = rmg_cell(M,R,K,T.cny);
    if(!X.ok || X.i>=T.cnx)
    return;
    const int I = X.i, J = X.j, k = K;

    const int i0 = T.rx*I, j0 = T.ry*J;
    const int i1 = (T.rx==2 && i0+1<T.fnx) ? i0+1 : -1;
    const int j1 = (T.ry==2 && j0+1<T.fny) ? j0+1 : -1;

    rmg_ijk y;
    y.k = k; y.ok = true;

    y.i = i0; y.j = j0; y.c = (i0+j0)&1; y.b = (j0+1)&1;
    const float r00 = rmg_resid_at(F,y,f,u,cp,cn,cs,cw,ce,ct,cb);
    float s = r00;

    if(i1>=0 && j1>=0)
    {
        y.i = i0; y.j = j1; y.c = (i0+j1)&1; y.b = (j1+1)&1;
        const float r01 = rmg_resid_at(F,y,f,u,cp,cn,cs,cw,ce,ct,cb);
        y.i = i1; y.j = j0; y.c = (i1+j0)&1; y.b = (j0+1)&1;
        const float r10 = rmg_resid_at(F,y,f,u,cp,cn,cs,cw,ce,ct,cb);
        y.i = i1; y.j = j1; y.c = (i1+j1)&1; y.b = (j1+1)&1;
        const float r11 = rmg_resid_at(F,y,f,u,cp,cn,cs,cw,ce,ct,cb);
        s = (((r00 + r01) + r10) + r11)*0.25f;
    }
    else if(i1>=0)
    {
        y.i = i1; y.j = j0; y.c = (i1+j0)&1; y.b = (j0+1)&1;
        s = (r00 + rmg_resid_at(F,y,f,u,cp,cn,cs,cw,ce,ct,cb))*0.5f;
    }
    else if(j1>=0)
    {
        y.i = i0; y.j = j1; y.c = (i0+j1)&1; y.b = (j1+1)&1;
        s = (r00 + rmg_resid_at(F,y,f,u,cp,cn,cs,cw,ce,ct,cb))*0.5f;
    }

    const int qc = rmg_q(I,J,k,T.chh,T.csi,T.csk);
    cf[qc] = s;
    cu[qc] = 0.0f;
}

//  Fine cell: u += cell-centred bilinear interpolation of the coarse u in
//  the plane, identity in z (reefmg_core::prolong_xy), masked by act.  The
//  coarse places: J and Ja = J -+ 1 are of the other parity, so the four
//  coarse values are at place (J+1)/2 and (Ja+1)/2 of their halves.
inline void rmg_prolong(const rmg_tp T, int m, int r, int k,
                        const device float *cu, device float *fu,
                        const device unsigned char *act)
{
    if(m>=T.fhh)
    return;
    const rmg_ijk x = rmg_cell(m,r,k,T.fny);
    if(!x.ok || x.i>=T.fnx)
    return;
    const int i = x.i, j = x.j;

    const int q = k*T.fsk + (i+1)*T.fsi + x.c*T.fhh + m;
    if(act[q]==0)
    return;

    int I = i, Ia = i, J = j, Ja = j;
    float wx = 1.0f, wy = 1.0f;

    if(T.rx==2)
    {
        I = i >> 1;
        Ia = (i&1) ? I+1 : I-1;
        wx = 0.75f;
    }
    if(T.ry==2)
    {
        J = j >> 1;
        Ja = (j&1) ? J+1 : J-1;
        wy = 0.75f;
    }
    Ia = Ia<0 ? 0 : (Ia>T.cnx-1 ? T.cnx-1 : Ia);
    Ja = Ja<0 ? 0 : (Ja>T.cny-1 ? T.cny-1 : Ja);

    const float w00 = wx*wy, w10 = (1.0f-wx)*wy, w01 = wx*(1.0f-wy), w11 = (1.0f-wx)*(1.0f-wy);

    const int bI  = k*T.csk + (I +1)*T.csi;
    const int bIa = k*T.csk + (Ia+1)*T.csi;
    const int pJ  = (J +1)>>1, pJa = (Ja+1)>>1;

    const float c00 = cu[bI  + ((I +J )&1)*T.chh + pJ ];
    const float c10 = cu[bIa + ((Ia+J )&1)*T.chh + pJ ];
    const float c01 = cu[bI  + ((I +Ja)&1)*T.chh + pJa];
    const float c11 = cu[bIa + ((Ia+Ja)&1)*T.chh + pJa];

    fu[q] += w00*c00 + w10*c10 + w01*c01 + w11*c11;
}

//  Coarse operator of one coarse cell, reefmg_core::coarsen_t in single
//  precision: the face coefficients recovered from the fine rows,
//  averaged over the fine faces of the coarse face and rebuilt with the
//  coarse widths; the row sum carried over.  Widths hx/hy have one ghost
//  entry each side, index i+1 for cell i.
inline void rmg_coarsen(const rmg_tp T, int M, int R, int K,
                        const device float *fp, const device float *fn, const device float *fs,
                        const device float *fw, const device float *fe, const device float *ft,
                        const device float *fb, const device unsigned char *fact,
                        const device float *fhx, const device float *fhy,
                        device float *cp, device float *cn, device float *cs,
                        device float *cw, device float *ce, device float *ct,
                        device float *cb, device unsigned char *cact,
                        const device float *chx, const device float *chy)
{
    if(M>=T.chh)
    return;
    const rmg_ijk x = rmg_cell(M,R,K,T.cny);
    if(!x.ok || x.i>=T.cnx)
    return;
    const int I = x.i, J = x.j, k = K;
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
        const int qf = rmg_q(i,j,k,T.fhh,T.fsi,T.fsk);
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

    const int qc = rmg_q(I,J,k,T.chh,T.csi,T.csk);

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
inline void rmg_factor(const rmg_lp P, int m, int r,
                       const device float *cp, const device float *ct, const device float *cb,
                       device float *ctc, device float *cti, const device unsigned char *act)
{
    if(m>=P.hh)
    return;
    const rmg_ijk x = rmg_cell(m,r,0,P.ny);
    if(!x.ok || x.i>=P.nx)
    return;
    const int sk = P.sk;
    const int q0 = rmg_q(x.i,x.j,0,P.hh,P.si,sk);

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

//  rmg_factor, and for rmg_smooth_h the couplings of every cell divided by
//  its pivot, ti = 1/pivot in single precision, hc = ctc
inline void rmg_factor_h(const rmg_lp P, int m, int r,
                         const device float *cp, const device float *cn, const device float *cs,
                         const device float *cw, const device float *ce, const device float *ct,
                         const device float *cb, device float *cti,
                         device rmg_h *an, device rmg_h *as, device rmg_h *aw, device rmg_h *ae,
                         device rmg_h *ab, device rmg_h *hc, const device unsigned char *act)
{
    if(m>=P.hh)
    return;
    const rmg_ijk x = rmg_cell(m,r,0,P.ny);
    if(!x.ok || x.i>=P.nx)
    return;
    const int sk = P.sk;
    const int q0 = rmg_q(x.i,x.j,0,P.hh,P.si,sk);

    bool any = false;
    for(int k=0; k<P.nz; ++k)
    if(act[q0+k*sk]!=0)
    {
        any = true;
        break;
    }
    if(!any)
    return;

    float tcp = 0.0f;
    for(int k=0; k<P.nz; ++k)
    {
        const int q = q0 + k*sk;
        const float inv = (k==0) ? 1.0f/cp[q] : 1.0f/(cp[q] - cb[q]*tcp);
        tcp = ct[q]*inv;
        cti[q] = inv;
        hc[q] = rmg_toh(tcp);
        an[q] = rmg_toh(cn[q]*inv);
        as[q] = rmg_toh(cs[q]*inv);
        aw[q] = rmg_toh(cw[q]*inv);
        ae[q] = rmg_toh(ce[q]*inv);
        ab[q] = rmg_toh(cb[q]*inv);
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
        u[rmg_q(-1 ,j,k,P.hh,P.si,P.sk)] = u[rmg_q(nx_-1,j,k,P.hh,P.si,P.sk)];
        u[rmg_q(nx_,j,k,P.hh,P.si,P.sk)] = u[rmg_q(0    ,j,k,P.hh,P.si,P.sk)];
        return;
    }

    const int iy = it - ny_*nz_;
    if(P.pery==0 || iy>=nx_*nz_)
    return;
    const int i = iy % nx_;
    const int k = iy / nx_;
    u[rmg_q(i,-1 ,k,P.hh,P.si,P.sk)] = u[rmg_q(i,ny_-1,k,P.hh,P.si,P.sk)];
    u[rmg_q(i,ny_,k,P.hh,P.si,P.sk)] = u[rmg_q(i,0    ,k,P.hh,P.si,P.sk)];
}

inline int rmg_nwrap(const rmg_lp P)
{
    return P.ny*P.nz + P.nx*P.nz;
}

inline void rmg_zero(unsigned t, device float *a)
{
    a[t] = 0.0f;
}

//  ------------------------------------------------------------ the tail
//  The small coarse levels lt .. nuse-1 of the hierarchy, one block of
//  memory from level lt on (B); every array of a level at a byte offset.
struct rmg_tlv
{
    rmg_lp P;
    rmg_tp T;     // to the next coarser level (not used by the last one)
    int of, ou, op, on, os, ow, oe, ot, ob, otc, oti, oact;
    int ohn, ohs, ohw, ohe, ohb, ohc;
};

struct rmg_tail
{
    int nlev, pre, post, dpre, dpost, csweeps, pad0, pad1;
    rmg_tlv L[12];
};

inline device float *rmg_pf(device char *B, int o)
{
    return (device float*)(B + o);
}

inline device rmg_h *rmg_ph(device char *B, int o)
{
    return (device rmg_h*)(B + o);
}

inline device unsigned char *rmg_pa(device char *B, int o)
{
    return (device unsigned char*)(B + o);
}

inline void rmg_tail_wrap(constant rmg_tlv &V, device char *B, unsigned tid, unsigned tg)
{
    const rmg_lp P = V.P;
    const int nw = rmg_nwrap(P);
    for(int t=int(tid); t<nw; t+=int(tg))
    rmg_wrap(P,unsigned(t),rmg_pf(B,V.ou));
    rmg_barrier();
}

//  reefmg_device::smooth: the periodic halo before every half-sweep
inline void rmg_tail_smooth(constant rmg_tlv &V, device char *B, int sweeps, int dir,
                            unsigned tid, unsigned tg)
{
    const rmg_lp P = V.P;
    const int n = P.nx*P.hh;
    const bool per = (P.perx!=0 || P.pery!=0);
    for(int sw=0; sw<sweeps; ++sw)
    for(int h=0; h<2; ++h)
    {
        if(per)
        rmg_tail_wrap(V,B,tid,tg);
        const int colour = (dir==1) ? 1-h : h;
        for(int t=int(tid); t<n; t+=int(tg))
        {
            const int m = t % P.hh, i = t / P.hh;
            if(P.h16!=0)
            rmg_smooth_h(P,colour,m,i,0,rmg_pf(B,V.of),rmg_pf(B,V.ou),
                         rmg_ph(B,V.ohn),rmg_ph(B,V.ohs),rmg_ph(B,V.ohw),rmg_ph(B,V.ohe),
                         rmg_ph(B,V.ohb),rmg_ph(B,V.ohc),rmg_pf(B,V.oti),rmg_pa(B,V.oact));
            else
            rmg_smooth(P,colour,m,i,0,rmg_pf(B,V.of),rmg_pf(B,V.ou),
                       rmg_pf(B,V.on),rmg_pf(B,V.os),rmg_pf(B,V.ow),rmg_pf(B,V.oe),
                       rmg_pf(B,V.ob),rmg_pf(B,V.otc),rmg_pf(B,V.oti),rmg_pa(B,V.oact));
        }
        rmg_barrier();
    }
}

//  reefmg_device::cycle from level lt (L[0]) down, without recursion (MSL
//  has none): work item tid of tg does every tg-th item of every step
inline void rmg_tail_run(constant rmg_tail &X, device char *B, unsigned tid, unsigned tg)
{
    const int nl = X.nlev;

    for(int l=0; l<nl-1; ++l)
    {
        rmg_tail_smooth(X.L[l],B,X.pre,X.dpre,tid,tg);
        const rmg_lp P = X.L[l].P;
        if(P.perx!=0 || P.pery!=0)
        rmg_tail_wrap(X.L[l],B,tid,tg);

        const rmg_tp T = X.L[l].T;
        const int nr = 2*T.cnx;
        const int n = T.chh*nr*T.nz;
        for(int t=int(tid); t<n; t+=int(tg))
        {
            const int M = t % T.chh, w = t / T.chh;
            rmg_restrict_residual(T,P,M,w%nr,w/nr,
                                  rmg_pf(B,X.L[l].of),rmg_pf(B,X.L[l].ou),rmg_pf(B,X.L[l].op),
                                  rmg_pf(B,X.L[l].on),rmg_pf(B,X.L[l].os),rmg_pf(B,X.L[l].ow),
                                  rmg_pf(B,X.L[l].oe),rmg_pf(B,X.L[l].ot),rmg_pf(B,X.L[l].ob),
                                  rmg_pf(B,X.L[l+1].of),rmg_pf(B,X.L[l+1].ou));
        }
        rmg_barrier();
    }

    rmg_tail_smooth(X.L[nl-1],B,X.csweeps,2,tid,tg);

    for(int l=nl-2; l>=0; --l)
    {
        const rmg_tp T = X.L[l].T;
        const int nr = 2*T.fnx;
        const int n = T.fhh*nr*T.nz;
        for(int t=int(tid); t<n; t+=int(tg))
        {
            const int m = t % T.fhh, w = t / T.fhh;
            rmg_prolong(T,m,w%nr,w/nr,rmg_pf(B,X.L[l+1].ou),rmg_pf(B,X.L[l].ou),rmg_pa(B,X.L[l].oact));
        }
        rmg_barrier();
        rmg_tail_smooth(X.L[l],B,X.post,X.dpost,tid,tg);
    }
}

)

#undef device
#undef constant

#endif
