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

//  The REEFMG preconditioner on a device (N 16).
//
//  A device holds ONE multigrid hierarchy of the WHOLE domain, in single
//  precision and in the device layout of reefmg_device_kernels.h, and runs
//  everything of the preconditioner on it: the coarse operators, the column
//  factorisations and the V-cycle.  reefmg_core keeps BiCGStab, the exact
//  operator and the convergence test on the CPU ranks, in double; the ranks
//  write their blocks of the right-hand side straight into the device's fine
//  level and read the correction back (an MPI shared-memory window on one
//  node, see reefmg_core::device_setup).
//
//  Work is queued in order and completed at sync(); the arrays returned by
//  arr() are host-visible (unified memory on Apple Silicon) and may be used
//  by the host only between sync() and the next queued kernel.
//
//  Devices:
//    1  Metal - the GPU of an Apple Silicon Mac (reefmg_device_metal.cpp,
//       in every macOS build: CodeLite, make or cmake, nothing to switch on)
//    2  reference - the same kernels in plain loops on the CPU, single
//       threaded; for testing the device path on any machine

//  arrays of a level
enum rmg_array {RMG_P=0, RMG_N, RMG_S, RMG_W, RMG_E, RMG_T, RMG_B, RMG_TC, RMG_TI,
                RMG_U, RMG_F, RMG_R, RMG_ACT, RMG_HX, RMG_HY, RMG_NARRAY};

//  Byte offsets of the arrays of one level in a single block of memory, each
//  aligned to 256 bytes; floats, except RMG_ACT (one byte per cell), with
//  the row layout of reefmg_device_kernels.h (the two colours of a row in
//  separate halves).  The host computes the same layout for the fine level
//  it shares with the device.
struct rmg_layout
{
    long off[RMG_NARRAY];
    long bytes;

    int hh,si;   // entries per colour of a row, row length (see reefmg_device_kernels.h)
    long sk;     // plane

    rmg_layout(int nx,int ny,int nz)
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
            o+=(n+255)/256*256;
        }
        bytes=o;
    }
};

struct rmg_level_desc
{
    int nx,ny,nz;   // extent
    int rx,ry;      // coarsening ratio from the level above (1 or 2), ignored for level 0
};

class reefmg_device
{
public:
    virtual ~reefmg_device(){}

    //  "Metal GPU Apple M3 Max", "reference (CPU, single threaded)"
    virtual std::string name() const=0;

    //  The hierarchy, level 0 first.  ext0: if not null, level 0 lives in
    //  this host memory (rmg_layout(level 0).bytes, aligned to the VM page
    //  size), shared with the caller; returns false with a message if the
    //  device cannot use it (the caller retries with ext0 = null and copies).
    virtual bool build(const std::vector<rmg_level_desc> &lv,int perx,int pery,void *ext0,std::string &err)=0;

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

    int levels() const {return (int)desc.size();}
    const rmg_level_desc& level(int l) const {return desc[l];}

protected:
    //  queued kernels, implemented by each device
    //  zero: first half-sweep from u = 0 (the other colour not read, set to 0)
    virtual void k_smooth(int l,int colour,int zero)=0;
    virtual void k_restrict(int l)=0;        // f(l+1) = R (f - A u)(l), u(l+1) = 0
    virtual void k_prolong(int l)=0;         // u(l) += P u(l+1)
    virtual void k_coarsen(int l)=0;         // level l+1 from level l
    virtual void k_factor(int l)=0;
    virtual void k_wrap(int l)=0;            // periodic halo of u(l)
    virtual void k_zero(int l,int a)=0;
    //  sweeps symmetric red-black sweeps on level l, as one unit of work
    //  where the device can (Metal: one threadgroup for a small level)
    virtual void k_coarse_solve(int l,int sweeps);

    void smooth(int l,int sweeps,int dir,bool from_zero=false);
    void cycle(int l,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,bool from_zero=false);

    std::vector<rmg_level_desc> desc;
    int perx=0, pery=0;
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
