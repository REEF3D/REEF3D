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

#ifndef REEFMG_DEVICE_H_
#define REEFMG_DEVICE_H_

#include <string>
#include <vector>

namespace rmg_kern {struct rmg_lp; struct rmg_tp; struct rmg_tail;}

//  The REEFMG preconditioner on a device (N 16).
//
//  A device holds ONE multigrid hierarchy of the WHOLE domain, in single
//  precision and in the device layout of reefmg_device_kernels.h, and runs
//  everything of the preconditioner on it: the coarse operators, the column
//  factorisations and the V-cycle - and with N 16 2 and 4 a whole Krylov
//  solve in single precision (krylov()).  reefmg_core keeps the outer
//  Krylov method (BiCGStab, or GCR around krylov()), the exact operator and
//  the convergence test on the CPU ranks, in double; the ranks
//  write their blocks of the right-hand side straight into the device's fine
//  level and read the correction back (an MPI shared-memory window on one
//  node, see reefmg_core::device_setup).
//
//  Work is queued in order and completed at sync(); the arrays returned by
//  arr() are host-visible (unified memory on Apple Silicon) and may be used
//  by the host only between sync() and the next queued kernel.
//
//  Devices (reefmg_device_create; N 16 1 and 2 the first, 3 and 4 the
//  second - 2 and 4 with the Krylov solve on the device, krylov()):
//    1  Metal - the GPU of an Apple Silicon Mac (reefmg_device_metal.cpp,
//       in every macOS build: CodeLite, make or cmake, nothing to switch on)
//    2  reference - the same kernels in plain loops on the CPU, single
//       threaded; for testing the device path on any machine

//  arrays of a level
//  RMG_HN .. RMG_HC: the smoother's couplings in half precision (N 21 16,
//  rmg_smooth_h), held instead of RMG_TC
//  RMG_KR .. RMG_KT: the vectors of the Krylov solve on the device (N 16 2
//  and 4, krylov()), level 0 only: residual, search direction, A y, the
//  preconditioned direction y, the preconditioned residual z, A z (the
//  shadow residual is RMG_F); RMG_KD: the partial sums of the dot products
//  (RMG_DOTG per product)
enum rmg_array {RMG_P=0, RMG_N, RMG_S, RMG_W, RMG_E, RMG_T, RMG_B, RMG_TC, RMG_TI,
                RMG_U, RMG_F, RMG_R, RMG_ACT, RMG_HX, RMG_HY,
                RMG_HN, RMG_HS, RMG_HW, RMG_HE, RMG_HB, RMG_HC,
                RMG_KR, RMG_KP, RMG_KV, RMG_KY, RMG_KZ, RMG_KT, RMG_KD, RMG_NARRAY};

//  work items of a dot product, each summing every RMG_DOTG-th entry
const int RMG_DOTG=16384;

//  Byte offsets of the arrays of one level in a single block of memory, each
//  aligned to 256 bytes; floats, except RMG_ACT (one byte per cell) and the
//  half-precision couplings (two bytes), with the row layout of
//  reefmg_device_kernels.h (the two colours of a row in separate halves).
//  h16: the smoother's coefficients in half precision (RMG_HN .. RMG_HC
//  instead of RMG_TC); kry: the vectors of the Krylov solve (RMG_KR ..
//  RMG_KD, level 0 of N 16 2 and 4).  The host computes the same layout for
//  the fine level it shares with the device.
struct rmg_layout
{
    long off[RMG_NARRAY];
    long bytes;

    int hh,si;   // entries per colour of a row, row length (see reefmg_device_kernels.h)
    long sk;     // plane

    rmg_layout(int nx,int ny,int nz,bool h16,bool kry=false)
    {
        hh=(ny+3)/2;
        si=2*hh;
        sk=(long)(nx+2)*si;
        const long N=sk*nz;
        long o=0;
        for(int a=0;a<RMG_NARRAY;++a)
        {
            off[a]=o;
            long n=N*4;
            if(a==RMG_ACT) n=N;
            if(a==RMG_R)   n=0;     // the residual is not stored (k_restrict computes it)
            if(a==RMG_HX)  n=(long)(nx+2)*4;
            if(a==RMG_HY)  n=(long)(ny+2)*4;
            if(a==RMG_TC && h16) n=0;
            if(a>=RMG_HN && a<=RMG_HC) n=h16? N*2 : 0;
            if(a>=RMG_KR && a<=RMG_KT) n=kry? N*4 : 0;
            if(a==RMG_KD) n=kry? 2L*RMG_DOTG*4 : 0;
            o+=(n+255)/256*256;
        }
        bytes=o;
    }
};

struct rmg_level_desc
{
    int nx,ny,nz;   // extent
    int rx,ry;      // coarsening ratio from the level above (1 or 2), ignored for level 0
    int kry;        // level 0: the vectors of the Krylov solve (krylov()), 0 or 1
};

class reefmg_device
{
public:
    reefmg_device();
    virtual ~reefmg_device(){}

    //  "Metal GPU Apple M3 Max", "reference (CPU, single threaded)"
    virtual std::string name() const=0;

    //  The hierarchy, level 0 first; h16: half-precision smoother
    //  coefficients (N 21).  ext0: if not null, level 0 lives in this host
    //  memory (rmg_layout(level 0).bytes, aligned to the VM page size),
    //  shared with the caller; returns false with a message if the device
    //  cannot use it (the caller retries with ext0 = null and copies).  The
    //  levels from 1 on are one block of memory (for the tail, below).
    virtual bool build(const std::vector<rmg_level_desc> &lv,int perx,int pery,bool h16,void *ext0,std::string &err)=0;

    //  host-visible array a of level l
    virtual void* arr(int l,int a)=0;

    //  coarse operators and the factorisations of levels 0 .. nuse-1 from
    //  the fine coefficients and activity (queued)
    void coarsen(int nuse);

    //  u(0) = V-cycle applied to f(0), from u(0) = 0 (queued).  Pre- and
    //  post-smoothing as reefmg_core with ordering 1 (red-black): red then
    //  black before the coarse correction, black then red after it; the
    //  coarsest used level gets coarse_sweeps symmetric sweeps.
    //  ncycle > 1: that many V-cycles, each starting from the result of the
    //  previous one (stationary multigrid on the fp32 operator)
    void vcycle(int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,int ncycle=1);

    //  The Krylov solve on the device (N 16 2 and 4): u(0) = A^-1 f(0) with
    //  the device's single-precision fine operator, by BiCGStab from u = 0,
    //  right-preconditioned by ncycle V-cycles (as vcycle), until the
    //  residual has fallen by tol or after maxit iterations.  The scalars
    //  are summed and kept in double on the host, which waits for the device
    //  at every dot product.  Needs the Krylov vectors on level 0
    //  (rmg_level_desc::kry).  Returns the iterations; relres the residual
    //  reduction reached (in single precision).  f(0) is not changed.
    int krylov(int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,int ncycle,
               double tol,int maxit,double &relres);

    //  wait for everything queued; false with a message if the device failed
    virtual bool sync(std::string &err)=0;

    //  bytes held by the device (without an external level 0)
    virtual long memory_bytes() const=0;

    //  seconds the device was busy (Metal: GPU start to end of every command
    //  buffer) since the last reset; 0 where not measured
    virtual double gpu_time() const {return 0.0;}
    virtual void reset_gpu_time(){}

    //  Metal with the environment variable REEFMG_DEVICE_PROFILE=1: every
    //  kernel is committed and timed on its own (slow, for analysis only);
    //  this reports the GPU time per kind of kernel
    virtual std::string profile_report() const {return "";}

    //  V-cycles run since the last reset (vcycle and krylov)
    long vcycles() const {return nvc;}
    void reset_vcycles(){nvc=0;}

    int levels() const {return (int)desc.size();}
    const rmg_level_desc& level(int l) const {return desc[l];}
    bool half() const {return h16;}

    //  The tail: the small coarse levels at the bottom of a V-cycle run as
    //  ONE unit of work (Metal: one dispatch of one threadgroup) instead of
    //  one dispatch per step - a dispatch costs a few microseconds however
    //  little it does, and the steps of a level of a few thousand columns do
    //  less than that.  It starts at the first level from 1 on whose columns
    //  of one colour, nx*hh, are at most tail_columns; the coarsest level
    //  alone if it has up to 4096.  Default 2048, environment variable
    //  REEFMG_DEVICE_TAIL to try others (0: the coarsest level only, -1: no
    //  tail, every step its own dispatch).
    void set_tail_columns(long c){tailcols=c;}
    long tail_columns() const {return tailcols;}
    //  first level of the tail of a hierarchy of nuse levels; nuse if none
    int tail_level(int nuse) const;

protected:
    //  queued kernels, implemented by each device
    //  zero: first half-sweep from u = 0 (the other colour not read, set to 0)
    virtual void k_smooth(int l,int colour,int zero)=0;
    virtual void k_restrict(int l)=0;        // f(l+1) = R (f - A u)(l), u(l+1) = 0
    virtual void k_prolong(int l)=0;         // u(l) += P u(l+1)
    virtual void k_coarsen(int l)=0;         // level l+1 from level l
    virtual void k_factor(int l)=0;
    virtual void k_wrap(int l,int a)=0;      // periodic halo of array a of level l
    virtual void k_zero(int l,int a)=0;
    //  the Krylov solve, level 0 (arrays of level l):
    //  y = A u (f < 0) or y = f - A u, at every cell of the level
    virtual void k_apply(int l,int f,int u,int y)=0;
    //  w = c0 w + c1 a + c2 b, every entry of the level (w not read if c0 = 0)
    virtual void k_axpy3(int l,int w,float c0,int a,float c1,int b,float c2)=0;
    //  the partial sums of (a,b) and (c,d) into RMG_KD
    virtual void k_dot2(int l,int a,int b,int c,int d)=0;
    //  k_axpy3, and the partial sums of (w,w) and (w,h) of the new w into RMG_KD
    virtual void k_axpy3_dot(int l,int w,float c0,int a,float c1,int b,float c2,int h)=0;
    //  the V-cycle of levels lt .. nuse-1 (lt >= 1) as one unit of work,
    //  from the parameters of tail_params; by default step by step
    virtual void k_tail(int lt,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps);

    //  the arrays the kernels use for f and u of level 0: RMG_F and RMG_U,
    //  or two of the Krylov vectors while krylov() applies the V-cycle to
    //  them; amap() gives the array a kernel binds for array a of level l
    int fa=RMG_F, ua=RMG_U;
    int amap(int l,int a) const
    {
        if(l==0 && a==RMG_F) return fa;
        if(l==0 && a==RMG_U) return ua;
        return a;
    }

    //  krylov(): y = V(x) with the V-cycle of vcycle(), and the two dot
    //  products queued last (sync, then the partial sums in double)
    void kprecond(int x,int y,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,int ncycle);
    bool kdots(double &d1,double &d2);

    void smooth(int l,int sweeps,int dir,bool from_zero=false);
    void cycle(int l,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,bool from_zero=false);

    //  build(): the layouts, and the byte offset of every level from 1 on
    //  in the block of the coarse levels; returns the block's size
    long plan_levels(const std::vector<rmg_level_desc> &lv,int px,int py,bool half16);
    //  the tail's parameters, offsets from the start of level lt
    void tail_params(int lt,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,rmg_kern::rmg_tail &X) const;
    rmg_kern::rmg_lp lp(int l) const;
    rmg_kern::rmg_tp tp(int l) const;

    std::vector<rmg_level_desc> desc;
    std::vector<rmg_layout> lay;
    std::vector<long> lbase;    // byte offset of level l in the coarse block (0 for level 0)
    int perx=0, pery=0;
    bool h16=false;
    long tailcols=2048;
    int tail0=-1;               // tail of the V-cycle being queued
    long nvc=0;                 // see vcycles()
};

//  kind 1 Metal, 2 reference.  Returns 0 with a message in err if the kind
//  is not available (not built in, no GPU).
reefmg_device* reefmg_device_create(int kind,std::string &err);

//  The Metal device (reefmg_device_metal.cpp) is part of every macOS build;
//  -DREEF3D_NO_METAL leaves it out.
#if defined(__APPLE__) && !defined(REEF3D_NO_METAL)
#define REEF3D_METAL 1
reefmg_device* reefmg_device_metal_create(std::string &err);
#else
#define REEF3D_METAL 0
#endif

#endif
