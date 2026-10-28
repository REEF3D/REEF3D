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

//  REEFMG device for the Apple Silicon GPU (N 16 1), Objective-C++ with ARC.
//
//  Built only on macOS with make METAL=1 / cmake -DREEF3D_USE_METAL=ON
//  (clang++ -x objective-c++ -fobjc-arc, -framework Metal -framework
//  Foundation), without the C++ precompiled header.
//
//  - The kernels are the text of reefmg_device_kernels.h plus the entry
//    points below, compiled at run time (newLibraryWithSource), once per
//    process and shared by every hierarchy: no offline Metal compiler.
//  - Every level is one MTLBuffer in shared storage, its arrays at the
//    offsets of rmg_layout.  The fine level can be the caller's memory (the
//    MPI shared window of reefmg_core), wrapped without a copy.
//  - Kernels are queued into one compute encoder of serial dispatch type and
//    committed at sync(): one V-cycle is one command buffer and one wait.
//    The coarsest level's sweeps run as one dispatch of one threadgroup.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "reefmg_device.h"
#include "reefmg_device_kernels.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>

namespace
{
    //  entry points; buffer(0) the parameters, buffer(1) a scalar, the
    //  arrays from buffer(2) on in the order of the bind lists below
    const char *rmg_metal_entries = R"MSL(

kernel void k_smooth(constant rmg_lp &P [[buffer(0)]], constant int &colour [[buffer(1)]],
                     const device float *f [[buffer(2)]], device float *u [[buffer(3)]],
                     const device float *cn [[buffer(4)]], const device float *cs [[buffer(5)]],
                     const device float *cw [[buffer(6)]], const device float *ce [[buffer(7)]],
                     const device float *cb [[buffer(8)]], const device float *ctc [[buffer(9)]],
                     const device float *cti [[buffer(10)]], const device unsigned char *act [[buffer(11)]],
                     uint t [[thread_position_in_grid]])
{
    rmg_smooth(P,colour,t,f,u,cn,cs,cw,ce,cb,ctc,cti,act);
}

kernel void k_coarse(constant rmg_lp &P [[buffer(0)]], constant int &sweeps [[buffer(1)]],
                     const device float *f [[buffer(2)]], device float *u [[buffer(3)]],
                     const device float *cn [[buffer(4)]], const device float *cs [[buffer(5)]],
                     const device float *cw [[buffer(6)]], const device float *ce [[buffer(7)]],
                     const device float *cb [[buffer(8)]], const device float *ctc [[buffer(9)]],
                     const device float *cti [[buffer(10)]], const device unsigned char *act [[buffer(11)]],
                     uint tid [[thread_position_in_threadgroup]], uint tg [[threads_per_threadgroup]])
{
    const uint nc = uint(rmg_ncolour(P));
    const uint nw = uint(rmg_nwrap(P));
    const bool per = (P.perx!=0 || P.pery!=0);

    for(int sw=0; sw<sweeps; ++sw)
    for(int h=0; h<2; ++h)
    {
        if(per)
        {
            for(uint t=tid; t<nw; t+=tg) rmg_wrap(P,t,u);
            threadgroup_barrier(mem_flags::mem_device);
        }
        for(uint t=tid; t<nc; t+=tg) rmg_smooth(P,h,t,f,u,cn,cs,cw,ce,cb,ctc,cti,act);
        threadgroup_barrier(mem_flags::mem_device);
    }
}

kernel void k_residual(constant rmg_lp &P [[buffer(0)]],
                       const device float *f [[buffer(2)]], const device float *u [[buffer(3)]],
                       device float *r [[buffer(4)]], const device float *cp [[buffer(5)]],
                       const device float *cn [[buffer(6)]], const device float *cs [[buffer(7)]],
                       const device float *cw [[buffer(8)]], const device float *ce [[buffer(9)]],
                       const device float *ct [[buffer(10)]], const device float *cb [[buffer(11)]],
                       uint t [[thread_position_in_grid]])
{
    rmg_residual(P,t,f,u,r,cp,cn,cs,cw,ce,ct,cb);
}

kernel void k_restrict(constant rmg_tp &T [[buffer(0)]],
                       const device float *fr [[buffer(2)]], device float *cf [[buffer(3)]],
                       device float *cu [[buffer(4)]],
                       uint t [[thread_position_in_grid]])
{
    rmg_restrict(T,t,fr,cf,cu);
}

kernel void k_prolong(constant rmg_tp &T [[buffer(0)]],
                      const device float *cu [[buffer(2)]], device float *fu [[buffer(3)]],
                      const device unsigned char *act [[buffer(4)]],
                      uint t [[thread_position_in_grid]])
{
    rmg_prolong(T,t,cu,fu,act);
}

kernel void k_coarsen(constant rmg_tp &T [[buffer(0)]],
                      const device float *fp [[buffer(2)]], const device float *fn [[buffer(3)]],
                      const device float *fs [[buffer(4)]], const device float *fw [[buffer(5)]],
                      const device float *fe [[buffer(6)]], const device float *ft [[buffer(7)]],
                      const device float *fb [[buffer(8)]], const device unsigned char *fact [[buffer(9)]],
                      const device float *fhx [[buffer(10)]], const device float *fhy [[buffer(11)]],
                      device float *cp [[buffer(12)]], device float *cn [[buffer(13)]],
                      device float *cs [[buffer(14)]], device float *cw [[buffer(15)]],
                      device float *ce [[buffer(16)]], device float *ct [[buffer(17)]],
                      device float *cb [[buffer(18)]], device unsigned char *cact [[buffer(19)]],
                      const device float *chx [[buffer(20)]], const device float *chy [[buffer(21)]],
                      uint t [[thread_position_in_grid]])
{
    rmg_coarsen(T,t,fp,fn,fs,fw,fe,ft,fb,fact,fhx,fhy,cp,cn,cs,cw,ce,ct,cb,cact,chx,chy);
}

kernel void k_factor(constant rmg_lp &P [[buffer(0)]],
                     const device float *cp [[buffer(2)]], const device float *ct [[buffer(3)]],
                     const device float *cb [[buffer(4)]], device float *ctc [[buffer(5)]],
                     device float *cti [[buffer(6)]], const device unsigned char *act [[buffer(7)]],
                     uint t [[thread_position_in_grid]])
{
    rmg_factor(P,t,cp,ct,cb,ctc,cti,act);
}

kernel void k_wrap(constant rmg_lp &P [[buffer(0)]], device float *u [[buffer(2)]],
                   uint t [[thread_position_in_grid]])
{
    rmg_wrap(P,t,u);
}

kernel void k_zero(device float *a [[buffer(2)]], uint t [[thread_position_in_grid]])
{
    rmg_zero(t,a);
}

)MSL";

    enum {PS_SMOOTH=0, PS_COARSE, PS_RESIDUAL, PS_RESTRICT, PS_PROLONG, PS_COARSEN, PS_FACTOR, PS_WRAP, PS_ZERO, PS_N};
    const char *ps_names[PS_N]={"k_smooth","k_coarse","k_residual","k_restrict","k_prolong",
                                "k_coarsen","k_factor","k_wrap","k_zero"};

    struct metal_ctx
    {
        id<MTLDevice> device;
        id<MTLCommandQueue> queue;
        id<MTLComputePipelineState> ps[PS_N];
        std::string err;
        bool ok;
    };

    //  device, queue and kernels: once per process
    metal_ctx& context()
    {
        static metal_ctx c;
        static bool init=false;
        if(init) return c;
        init=true;

        @autoreleasepool
        {
            c.ok=true;
            c.device=MTLCreateSystemDefaultDevice();
            if(c.device==nil)
            {
                c.err="no Metal GPU found";
                c.ok=false;
                return c;
            }
            c.queue=[c.device newCommandQueue];

            std::string src="#include <metal_stdlib>\nusing namespace metal;\n";
            src+=rmg_kernel_source;
            src+="\n";
            src+=rmg_metal_entries;

            NSError *e=nil;
            id<MTLLibrary> lib=[c.device newLibraryWithSource:[NSString stringWithUTF8String:src.c_str()]
                                                      options:[MTLCompileOptions new]
                                                        error:&e];
            if(lib==nil)
            {
                c.err=std::string("Metal kernel compilation failed: ")
                     +(e? [[e localizedDescription] UTF8String] : "unknown error");
                c.ok=false;
                return c;
            }

            for(int k=0;k<PS_N && c.ok;++k)
            {
                id<MTLFunction> fn=[lib newFunctionWithName:[NSString stringWithUTF8String:ps_names[k]]];
                if(fn==nil)
                {
                    c.err=std::string("Metal kernel ")+ps_names[k]+" not found";
                    c.ok=false;
                    break;
                }
                NSError *pe=nil;
                c.ps[k]=[c.device newComputePipelineStateWithFunction:fn error:&pe];
                if(c.ps[k]==nil)
                {
                    c.err=std::string("Metal pipeline ")+ps_names[k]+": "
                         +(pe? [[pe localizedDescription] UTF8String] : "unknown error");
                    c.ok=false;
                }
            }
        }
        return c;
    }

    class reefmg_device_metal final : public reefmg_device
    {
    public:

        explicit reefmg_device_metal(metal_ctx &c) : ctx(c), cb(nil), enc(nil), bytes(0), failed(false)
        {
            const char *e=getenv("REEFMG_DEVICE_PROFILE");
            profile=(e && e[0]=='1');
            for(int k=0;k<PS_N;++k){prof[k]=0.0; profn[k]=0;}
        }

        ~reefmg_device_metal() override
        {
            std::string e;
            sync(e);
        }

        std::string name() const override
        {
            return std::string("Metal GPU ")+[[ctx.device name] UTF8String];
        }

        bool build(const std::vector<rmg_level_desc> &lv,int px,int py,void *ext0,std::string &err) override
        {
            desc=lv; perx=px; pery=py;
            buf.clear();
            lay.clear();
            bytes=0;

            for(size_t l=0;l<lv.size();++l)
            {
                if((long)(lv[l].nx+2)*(lv[l].ny+2)*lv[l].nz>=2147483647L)
                {
                    err="level too large for 32-bit indexing";
                    return false;
                }
                lay.push_back(rmg_layout(lv[l].nx,lv[l].ny,lv[l].nz));
                const long len=lay.back().bytes;
                id<MTLBuffer> b=nil;

                if(l==0 && ext0)
                {
                    const long page=sysconf(_SC_PAGESIZE)>0? sysconf(_SC_PAGESIZE) : 16384;
                    const long wlen=(len+page-1)/page*page;
                    b=[ctx.device newBufferWithBytesNoCopy:ext0 length:(NSUInteger)wlen
                                                   options:MTLResourceStorageModeShared deallocator:nil];
                    if(b==nil)
                    {
                        err="Metal cannot use the shared fine level without a copy";
                        return false;
                    }
                }
                else
                {
                    b=[ctx.device newBufferWithLength:(NSUInteger)len options:MTLResourceStorageModeShared];
                    if(b==nil)
                    {
                        err="Metal buffer allocation of "+std::to_string(len)+" bytes failed";
                        return false;
                    }
                    std::memset([b contents],0,(size_t)len);
                    bytes+=len;
                }
                buf.push_back(b);
            }
            return true;
        }

        void* arr(int l,int a) override
        {
            return (char*)[buf[l] contents]+lay[l].off[a];
        }

        bool sync(std::string &err) override
        {
            if(failed)
            {
                err=fail;
                return false;
            }
            if(enc==nil) return true;

            @autoreleasepool
            {
                [enc endEncoding];
                [cb commit];
                [cb waitUntilCompleted];
                gpu+=[cb GPUEndTime]-[cb GPUStartTime];

                if([cb status]==MTLCommandBufferStatusError)
                {
                    NSError *e=[cb error];
                    fail=std::string("command buffer failed: ")+(e? [[e localizedDescription] UTF8String] : "unknown error");
                    failed=true;
                }
                enc=nil;
                cb=nil;
            }

            if(failed)
            {
                err=fail;
                return false;
            }
            return true;
        }

        long memory_bytes() const override {return bytes;}

        double gpu_time() const override {return gpu;}
        void reset_gpu_time() override {gpu=0.0;}

        std::string profile_report() const override
        {
            if(!profile) return "";
            std::string r;
            double tot=0.0;
            for(int k=0;k<PS_N;++k) tot+=prof[k];
            char b[160];
            for(int k=0;k<PS_N;++k)
            if(profn[k]>0)
            {
                snprintf(b,sizeof(b),"    %-11s %8ld launches %9.3f ms  %5.1f %%\n",ps_names[k],profn[k],1.0e3*prof[k],
                         tot>0.0? 100.0*prof[k]/tot : 0.0);
                r+=b;
            }
            return r;
        }

    protected:

        void k_smooth(int l,int colour) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const long n=rmg_kern::rmg_ncolour(P);
            if(n<=0) return;
            id<MTLComputeCommandEncoder> e=begin(PS_SMOOTH);
            [e setBytes:&P length:sizeof(P) atIndex:0];
            [e setBytes:&colour length:sizeof(int) atIndex:1];
            const int a[]={RMG_F,RMG_U,RMG_N,RMG_S,RMG_W,RMG_E,RMG_B,RMG_TC,RMG_TI,RMG_ACT};
            bind(e,l,a,10,2);
            dispatch(e,PS_SMOOTH,n);
        }

        //  the sweeps of a small coarsest level in one threadgroup
        void k_coarse_solve(int l,int sweeps) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const long nc=rmg_kern::rmg_ncolour(P);
            if(nc>4096)
            {
                reefmg_device::k_coarse_solve(l,sweeps);
                return;
            }
            id<MTLComputeCommandEncoder> e=begin(PS_COARSE);
            [e setBytes:&P length:sizeof(P) atIndex:0];
            [e setBytes:&sweeps length:sizeof(int) atIndex:1];
            const int a[]={RMG_F,RMG_U,RMG_N,RMG_S,RMG_W,RMG_E,RMG_B,RMG_TC,RMG_TI,RMG_ACT};
            bind(e,l,a,10,2);
            NSUInteger w=[ctx.ps[PS_COARSE] maxTotalThreadsPerThreadgroup];
            if(w>1024) w=1024;
            [e dispatchThreadgroups:MTLSizeMake(1,1,1) threadsPerThreadgroup:MTLSizeMake(w,1,1)];
            if(profile) flush_profile(PS_COARSE);
        }

        void k_residual(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            id<MTLComputeCommandEncoder> e=begin(PS_RESIDUAL);
            [e setBytes:&P length:sizeof(P) atIndex:0];
            const int a[]={RMG_F,RMG_U,RMG_R,RMG_P,RMG_N,RMG_S,RMG_W,RMG_E,RMG_T,RMG_B};
            bind(e,l,a,10,2);
            dispatch(e,PS_RESIDUAL,(long)P.nx*P.ny*P.nz);
        }

        void k_restrict(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            id<MTLComputeCommandEncoder> e=begin(PS_RESTRICT);
            [e setBytes:&T length:sizeof(T) atIndex:0];
            const int af[]={RMG_R};
            const int ac[]={RMG_F,RMG_U};
            bind(e,l,af,1,2);
            bind(e,l+1,ac,2,3);
            dispatch(e,PS_RESTRICT,(long)T.cnx*T.cny*T.nz);
        }

        void k_prolong(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            id<MTLComputeCommandEncoder> e=begin(PS_PROLONG);
            [e setBytes:&T length:sizeof(T) atIndex:0];
            const int ac[]={RMG_U};
            const int af[]={RMG_U,RMG_ACT};
            bind(e,l+1,ac,1,2);
            bind(e,l,af,2,3);
            dispatch(e,PS_PROLONG,(long)T.fnx*T.fny*T.nz);
        }

        void k_coarsen(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            id<MTLComputeCommandEncoder> e=begin(PS_COARSEN);
            [e setBytes:&T length:sizeof(T) atIndex:0];
            const int af[]={RMG_P,RMG_N,RMG_S,RMG_W,RMG_E,RMG_T,RMG_B,RMG_ACT,RMG_HX,RMG_HY};
            bind(e,l,af,10,2);
            bind(e,l+1,af,10,12);
            dispatch(e,PS_COARSEN,(long)T.cnx*T.cny*T.nz);
        }

        void k_factor(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            id<MTLComputeCommandEncoder> e=begin(PS_FACTOR);
            [e setBytes:&P length:sizeof(P) atIndex:0];
            const int a[]={RMG_P,RMG_T,RMG_B,RMG_TC,RMG_TI,RMG_ACT};
            bind(e,l,a,6,2);
            dispatch(e,PS_FACTOR,(long)P.nx*P.ny);
        }

        void k_wrap(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            id<MTLComputeCommandEncoder> e=begin(PS_WRAP);
            [e setBytes:&P length:sizeof(P) atIndex:0];
            const int a[]={RMG_U};
            bind(e,l,a,1,2);
            dispatch(e,PS_WRAP,rmg_kern::rmg_nwrap(P));
        }

        void k_zero(int l,int v) override
        {
            id<MTLComputeCommandEncoder> e=begin(PS_ZERO);
            const int a[]={v};
            bind(e,l,a,1,2);
            dispatch(e,PS_ZERO,(long)(desc[l].nx+2)*(desc[l].ny+2)*desc[l].nz);
        }

    private:

        id<MTLComputeCommandEncoder> begin(int k)
        {
            if(enc==nil)
            {
                @autoreleasepool
                {
                    cb=[ctx.queue commandBuffer];
                    enc=[cb computeCommandEncoder];
                }
            }
            [enc setComputePipelineState:ctx.ps[k]];
            return enc;
        }

        //  arrays a[0..n-1] of level l at buffer indices first, first+1, ...
        void bind(id<MTLComputeCommandEncoder> e,int l,const int *a,int n,int first)
        {
            //  buf[l] holds the reference; the method takes an unretained array
            __unsafe_unretained id<MTLBuffer> b[16];
            NSUInteger o[16];
            for(int m=0;m<n;++m){b[m]=buf[l]; o[m]=(NSUInteger)lay[l].off[a[m]];}
            [e setBuffers:b offsets:o withRange:NSMakeRange((NSUInteger)first,(NSUInteger)n)];
        }

        void dispatch(id<MTLComputeCommandEncoder> e,int k,long n)
        {
            if(n<=0) return;
            NSUInteger w=[ctx.ps[k] maxTotalThreadsPerThreadgroup];
            if(w>256) w=256;
            [e dispatchThreads:MTLSizeMake((NSUInteger)n,1,1) threadsPerThreadgroup:MTLSizeMake(w,1,1)];
            if(profile) flush_profile(k);
        }

        //  profiling: this kernel alone in a command buffer, timed on the GPU
        void flush_profile(int k)
        {
            @autoreleasepool
            {
                [enc endEncoding];
                [cb commit];
                [cb waitUntilCompleted];
                const double t=[cb GPUEndTime]-[cb GPUStartTime];
                prof[k]+=t;
                ++profn[k];
                gpu+=t;
                enc=nil;
                cb=nil;
            }
        }

        rmg_kern::rmg_lp lp(int l) const
        {
            const rmg_level_desc &d=desc[l];
            rmg_kern::rmg_lp P;
            P.nx=d.nx; P.ny=d.ny; P.nz=d.nz; P.si=d.ny+2; P.sk=(d.nx+2)*(d.ny+2);
            P.perx=perx; P.pery=pery;
            return P;
        }

        rmg_kern::rmg_tp tp(int l) const
        {
            const rmg_level_desc &f=desc[l], &c=desc[l+1];
            rmg_kern::rmg_tp T;
            T.fnx=f.nx; T.fny=f.ny; T.nz=f.nz; T.fsi=f.ny+2; T.fsk=(f.nx+2)*(f.ny+2);
            T.cnx=c.nx; T.cny=c.ny; T.csi=c.ny+2; T.csk=(c.nx+2)*(c.ny+2);
            T.rx=c.rx; T.ry=c.ry; T.perx=perx; T.pery=pery;
            return T;
        }

        metal_ctx &ctx;
        std::vector<id<MTLBuffer> > buf;
        std::vector<rmg_layout> lay;
        id<MTLCommandBuffer> cb;
        id<MTLComputeCommandEncoder> enc;
        long bytes;
        bool failed;
        std::string fail;
        double gpu=0.0;
        bool profile;
        double prof[PS_N];
        long profn[PS_N];
    };
}

reefmg_device* reefmg_device_metal_create(std::string &err)
{
    metal_ctx &c=context();
    if(!c.ok)
    {
        err=c.err;
        return 0;
    }
    return new reefmg_device_metal(c);
}
