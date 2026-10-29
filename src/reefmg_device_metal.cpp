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

//  REEFMG device for the Apple Silicon GPU (N 16 1), in plain C++.
//
//  Compiled into every build; active on macOS (reefmg_device.h: REEF3D_METAL,
//  off with -DREEF3D_NO_METAL), elsewhere this file is empty.  No build
//  setting, no framework to link and no Objective-C compiler is needed - the
//  CodeLite project, the Makefile and CMake build it like any other source,
//  precompiled header included:
//
//  - The Objective-C runtime (libobjc) and the Metal framework are opened
//    with dlopen() when N 16 1 asks for the GPU, and Metal is driven through
//    objc_msgSend with the method names of its Objective-C API (the way
//    Apple's metal-cpp does it).  Objects are reference counted by hand:
//    what new.../alloc returns is released, what a command buffer hands out
//    is retained while it is used.
//  - The kernels are the text of reefmg_device_kernels.h plus the entry
//    points below, compiled at run time (newLibraryWithSource), once per
//    process and shared by every hierarchy: no offline Metal compiler.
//  - The fine level is one MTLBuffer in shared storage, its arrays at the
//    offsets of rmg_layout; it can be the caller's memory (the MPI shared
//    window of reefmg_core), wrapped without a copy.  The coarse levels are
//    one more buffer, every level at its offset in it.
//  - Kernels are queued into one compute encoder of serial dispatch type and
//    committed at sync(): one V-cycle is one command buffer and one wait.
//    The small coarse levels at the bottom of the V-cycle run as one
//    dispatch of one threadgroup (k_tail, see reefmg_device::tail_level).

#include "reefmg_device.h"

#if REEF3D_METAL

#include "reefmg_device_kernels.h"

#include <dlfcn.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace
{
    //  the four names reefmg_device_kernels.h leaves to the MSL side
    const char *rmg_metal_prelude = R"MSL(

typedef half rmg_h;
inline float rmg_hf(rmg_h x){return float(x);}
inline rmg_h rmg_toh(float x){return rmg_h(x);}
inline void rmg_barrier(){threadgroup_barrier(mem_flags::mem_device);}

)MSL";

    const char *rmg_metal_entries = R"MSL(

kernel void k_smooth(constant rmg_lp &P [[buffer(0)]], constant int2 &cz [[buffer(1)]],
                     const device float *f [[buffer(2)]], device float *u [[buffer(3)]],
                     const device float *cn [[buffer(4)]], const device float *cs [[buffer(5)]],
                     const device float *cw [[buffer(6)]], const device float *ce [[buffer(7)]],
                     const device float *cb [[buffer(8)]], const device float *ctc [[buffer(9)]],
                     const device float *cti [[buffer(10)]], const device unsigned char *act [[buffer(11)]],
                     uint2 g [[thread_position_in_grid]])
{
    rmg_smooth(P,cz.x,int(g.x),int(g.y),cz.y,f,u,cn,cs,cw,ce,cb,ctc,cti,act);
}

kernel void k_smooth_h(constant rmg_lp &P [[buffer(0)]], constant int2 &cz [[buffer(1)]],
                       const device float *f [[buffer(2)]], device float *u [[buffer(3)]],
                       const device half *an [[buffer(4)]], const device half *as [[buffer(5)]],
                       const device half *aw [[buffer(6)]], const device half *ae [[buffer(7)]],
                       const device half *ab [[buffer(8)]], const device half *hc [[buffer(9)]],
                       const device float *cti [[buffer(10)]], const device unsigned char *act [[buffer(11)]],
                       uint2 g [[thread_position_in_grid]])
{
    rmg_smooth_h(P,cz.x,int(g.x),int(g.y),cz.y,f,u,an,as,aw,ae,ab,hc,cti,act);
}

kernel void k_tail(constant rmg_tail &X [[buffer(0)]], device char *B [[buffer(1)]],
                   uint tid [[thread_position_in_threadgroup]], uint tg [[threads_per_threadgroup]])
{
    rmg_tail_run(X,B,tid,tg);
}

kernel void k_restrict(constant rmg_tp &T [[buffer(0)]], constant rmg_lp &F [[buffer(1)]],
                       const device float *f [[buffer(2)]], const device float *u [[buffer(3)]],
                       const device float *cp [[buffer(4)]], const device float *cn [[buffer(5)]],
                       const device float *cs [[buffer(6)]], const device float *cw [[buffer(7)]],
                       const device float *ce [[buffer(8)]], const device float *ct [[buffer(9)]],
                       const device float *cb [[buffer(10)]],
                       device float *cf [[buffer(11)]], device float *cu [[buffer(12)]],
                       uint3 g [[thread_position_in_grid]])
{
    rmg_restrict_residual(T,F,int(g.x),int(g.y),int(g.z),f,u,cp,cn,cs,cw,ce,ct,cb,cf,cu);
}

kernel void k_prolong(constant rmg_tp &T [[buffer(0)]],
                      const device float *cu [[buffer(2)]], device float *fu [[buffer(3)]],
                      const device unsigned char *act [[buffer(4)]],
                      uint3 g [[thread_position_in_grid]])
{
    rmg_prolong(T,int(g.x),int(g.y),int(g.z),cu,fu,act);
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
                      uint3 g [[thread_position_in_grid]])
{
    rmg_coarsen(T,int(g.x),int(g.y),int(g.z),fp,fn,fs,fw,fe,ft,fb,fact,fhx,fhy,cp,cn,cs,cw,ce,ct,cb,cact,chx,chy);
}

kernel void k_factor(constant rmg_lp &P [[buffer(0)]],
                     const device float *cp [[buffer(2)]], const device float *ct [[buffer(3)]],
                     const device float *cb [[buffer(4)]], device float *ctc [[buffer(5)]],
                     device float *cti [[buffer(6)]], const device unsigned char *act [[buffer(7)]],
                     uint2 g [[thread_position_in_grid]])
{
    rmg_factor(P,int(g.x),int(g.y),cp,ct,cb,ctc,cti,act);
}

kernel void k_factor_h(constant rmg_lp &P [[buffer(0)]],
                       const device float *cp [[buffer(2)]], const device float *cn [[buffer(3)]],
                       const device float *cs [[buffer(4)]], const device float *cw [[buffer(5)]],
                       const device float *ce [[buffer(6)]], const device float *ct [[buffer(7)]],
                       const device float *cb [[buffer(8)]], device float *cti [[buffer(9)]],
                       device half *an [[buffer(10)]], device half *as [[buffer(11)]],
                       device half *aw [[buffer(12)]], device half *ae [[buffer(13)]],
                       device half *ab [[buffer(14)]], device half *hc [[buffer(15)]],
                       const device unsigned char *act [[buffer(16)]],
                       uint2 g [[thread_position_in_grid]])
{
    rmg_factor_h(P,int(g.x),int(g.y),cp,cn,cs,cw,ce,ct,cb,cti,an,as,aw,ae,ab,hc,act);
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

    enum {PS_SMOOTH=0, PS_SMOOTH_H, PS_TAIL, PS_RESTRICT, PS_PROLONG, PS_COARSEN, PS_FACTOR, PS_FACTOR_H,
          PS_WRAP, PS_ZERO, PS_N};
    const char *ps_names[PS_N]={"k_smooth","k_smooth_h","k_tail","k_restrict","k_prolong",
                                "k_coarsen","k_factor","k_factor_h","k_wrap","k_zero"};

    //  ------------------------------------------- Objective-C runtime, by hand
    typedef void* obj;                 // id
    typedef void* sel;                 // SEL
    typedef unsigned long nsuint;      // NSUInteger
    struct mtl_size {nsuint w, h, d;}; // MTLSize
    struct ns_range {nsuint loc, len;};// NSRange

    const nsuint MTL_STORAGE_SHARED=0;     // MTLResourceStorageModeShared
    const long   MTL_CB_STATUS_ERROR=5;    // MTLCommandBufferStatusError

    struct objc_rt
    {
        bool ok=false;
        std::string err;
        obj (*get_class)(const char*)=nullptr;
        sel (*sel_reg)(const char*)=nullptr;
        void (*msg_send)()=nullptr;
        void* (*pool_push)()=nullptr;
        void (*pool_pop)(void*)=nullptr;
        obj (*create_device)()=nullptr;
    };

    objc_rt& rt()
    {
        static objc_rt r;
        static bool init=false;
        if(init) return r;
        init=true;

        void *ho=dlopen("/usr/lib/libobjc.A.dylib",RTLD_NOW|RTLD_GLOBAL);
        void *hf=dlopen("/System/Library/Frameworks/Foundation.framework/Foundation",RTLD_NOW|RTLD_GLOBAL);
        void *hm=dlopen("/System/Library/Frameworks/Metal.framework/Metal",RTLD_NOW|RTLD_GLOBAL);
        if(!ho || !hf || !hm)
        {
            r.err=std::string("cannot load ")+(!ho? "libobjc" : (!hf? "Foundation" : "Metal"))+": "+(dlerror()? dlerror() : "");
            return r;
        }
        r.get_class=(obj(*)(const char*))dlsym(ho,"objc_getClass");
        r.sel_reg=(sel(*)(const char*))dlsym(ho,"sel_registerName");
        r.msg_send=(void(*)())dlsym(ho,"objc_msgSend");
        r.pool_push=(void*(*)())dlsym(ho,"objc_autoreleasePoolPush");
        r.pool_pop=(void(*)(void*))dlsym(ho,"objc_autoreleasePoolPop");
        r.create_device=(obj(*)())dlsym(hm,"MTLCreateSystemDefaultDevice");
        r.ok=(r.get_class && r.sel_reg && r.msg_send && r.pool_push && r.pool_pop && r.create_device);
        if(!r.ok) r.err="Objective-C runtime or Metal entry points not found";
        return r;
    }

    //  [o name:args...] - objc_msgSend called through a pointer of the exact
    //  type of the method, as the arm64 calling convention requires
    template<class R,class... A>
    R msg(obj o,sel s,A... a)
    {
        typedef R (*fn)(obj,sel,A...);
        return ((fn)rt().msg_send)(o,s,a...);
    }

    //  selector, registered once per call site
    #define RMG_SEL(name) ([]{ static sel s_=rt().sel_reg(name); return s_; }())

    void release(obj o){if(o) msg<void>(o,RMG_SEL("release"));}
    obj  retain(obj o){return o? msg<obj>(o,RMG_SEL("retain")) : o;}

    //  autoreleased NSString
    obj nsstring(const char *s)
    {
        return msg<obj>(rt().get_class("NSString"),RMG_SEL("stringWithUTF8String:"),s);
    }

    std::string cstr(obj nsstr)
    {
        if(!nsstr) return std::string();
        const char *c=msg<const char*>(nsstr,RMG_SEL("UTF8String"));
        return c? std::string(c) : std::string();
    }

    std::string describe(obj nserror)
    {
        return nserror? cstr(msg<obj>(nserror,RMG_SEL("localizedDescription"))) : std::string("unknown error");
    }

    struct pool
    {
        void *p;
        pool() : p(rt().pool_push()) {}
        ~pool(){rt().pool_pop(p);}
    };

    //  ------------------------------------------------------------ context
    struct metal_ctx
    {
        obj device=nullptr;
        obj queue=nullptr;
        obj ps[PS_N]={};
        std::string err;
        bool ok=false;
    };

    //  device, queue and kernels: once per process
    metal_ctx& context()
    {
        static metal_ctx c;
        static bool init=false;
        if(init) return c;
        init=true;

        objc_rt &r=rt();
        if(!r.ok)
        {
            c.err=r.err;
            return c;
        }

        pool p;
        c.device=r.create_device();
        if(!c.device)
        {
            c.err="no Metal GPU found";
            return c;
        }
        c.queue=msg<obj>(c.device,RMG_SEL("newCommandQueue"));

        std::string src="#include <metal_stdlib>\nusing namespace metal;\n";
        src+=rmg_metal_prelude;
        src+="\n";
        src+=rmg_kernel_source;
        src+="\n";
        src+=rmg_metal_entries;

        obj opts=msg<obj>(msg<obj>(r.get_class("MTLCompileOptions"),RMG_SEL("alloc")),RMG_SEL("init"));
        obj e=nullptr;
        obj lib=msg<obj>(c.device,RMG_SEL("newLibraryWithSource:options:error:"),nsstring(src.c_str()),opts,&e);
        release(opts);
        if(!lib)
        {
            c.err="Metal kernel compilation failed: "+describe(e);
            return c;
        }

        c.ok=true;
        for(int k=0;k<PS_N && c.ok;++k)
        {
            obj fn=msg<obj>(lib,RMG_SEL("newFunctionWithName:"),nsstring(ps_names[k]));
            if(!fn)
            {
                c.err=std::string("Metal kernel ")+ps_names[k]+" not found";
                c.ok=false;
                break;
            }
            obj pe=nullptr;
            c.ps[k]=msg<obj>(c.device,RMG_SEL("newComputePipelineStateWithFunction:error:"),fn,&pe);
            release(fn);
            if(!c.ps[k])
            {
                c.err=std::string("Metal pipeline ")+ps_names[k]+": "+describe(pe);
                c.ok=false;
            }
        }
        release(lib);
        return c;
    }

    nsuint max_threads(obj pipeline)
    {
        return msg<nsuint>(pipeline,RMG_SEL("maxTotalThreadsPerThreadgroup"));
    }

    //  -------------------------------------------------------------- device
    class reefmg_device_metal final : public reefmg_device
    {
    public:

        explicit reefmg_device_metal(metal_ctx &c) : ctx(c)
        {
            const char *e=getenv("REEFMG_DEVICE_PROFILE");
            profile=(e && e[0]=='1');
            for(int k=0;k<PS_N;++k){prof[k]=0.0; profn[k]=0;}
        }

        ~reefmg_device_metal() override
        {
            std::string e;
            sync(e);
            free_buffers();
        }

        std::string name() const override
        {
            pool p;
            return "Metal GPU "+cstr(msg<obj>(ctx.device,RMG_SEL("name")));
        }

        bool build(const std::vector<rmg_level_desc> &lv,int px,int py,bool half16,void *ext0,std::string &err) override
        {
            free_buffers();
            const long nb=plan_levels(lv,px,py,half16);
            bytes=0;
            for(size_t l=0;l<lv.size();++l)
            if(lay[l].sk*lv[l].nz>=2147483647L)
            {
                err="level too large for 32-bit indexing";
                return false;
            }

            obj b=nullptr;
            const long len=lay[0].bytes;
            if(ext0)
            {
                const long page=sysconf(_SC_PAGESIZE)>0? sysconf(_SC_PAGESIZE) : 16384;
                const long wlen=(len+page-1)/page*page;
                b=msg<obj>(ctx.device,RMG_SEL("newBufferWithBytesNoCopy:length:options:deallocator:"),
                           ext0,(nsuint)wlen,MTL_STORAGE_SHARED,(void*)nullptr);
                if(!b)
                {
                    err="Metal cannot use the shared fine level without a copy";
                    return false;
                }
            }
            else
            {
                b=new_buffer(len,err);
                if(!b) return false;
                bytes+=len;
            }
            buf.assign(lv.size(),b);

            if(nb>0)
            {
                blk=new_buffer(nb,err);
                if(!blk) return false;
                bytes+=nb;
                for(size_t l=1;l<lv.size();++l) buf[l]=blk;
            }
            profl.assign(lv.size(),0.0);
            return true;
        }

        void* arr(int l,int a) override
        {
            return (char*)contents(buf[l])+lbase[l]+lay[l].off[a];
        }

        bool sync(std::string &err) override
        {
            if(failed)
            {
                err=fail;
                return false;
            }
            if(!enc) return true;

            commit_and_wait();

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
            //  by level; the tail is counted with its first level
            for(size_t l=0;l<profl.size();++l)
            if(profl[l]>0.0)
            {
                snprintf(b,sizeof(b),"    level %-2d %5d x %5d x %3d  %9.3f ms  %5.1f %%%s\n",(int)l,desc[l].nx,desc[l].ny,
                         desc[l].nz,1.0e3*profl[l],tot>0.0? 100.0*profl[l]/tot : 0.0,
                         ((int)l==proftail)? "  (tail: this level and below)" : "");
                r+=b;
            }
            return r;
        }

    protected:

        void k_smooth(int l,int colour,int zero) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int ps=h16? PS_SMOOTH_H : PS_SMOOTH;
            obj e=begin(ps,l);
            const int cz[2]={colour,zero};
            bytes_at(e,&P,sizeof(P),0);
            bytes_at(e,cz,sizeof(cz),1);
            const int a32[]={RMG_F,RMG_U,RMG_N,RMG_S,RMG_W,RMG_E,RMG_B,RMG_TC,RMG_TI,RMG_ACT};
            const int a16[]={RMG_F,RMG_U,RMG_HN,RMG_HS,RMG_HW,RMG_HE,RMG_HB,RMG_HC,RMG_TI,RMG_ACT};
            bind(e,l,h16? a16 : a32,10,2);
            dispatch3(e,ps,P.hh,P.nx,1);
        }

        //  levels lt .. nuse-1 in one threadgroup, from the block of the
        //  coarse levels at level lt
        void k_tail(int lt,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps) override
        {
            rmg_kern::rmg_tail X;
            tail_params(lt,nuse,pre,post,sweepstyle,coarse_sweeps,X);
            obj e=begin(PS_TAIL,lt);
            bytes_at(e,&X,sizeof(X),0);
            msg<void>(e,RMG_SEL("setBuffer:offset:atIndex:"),blk,(nsuint)lbase[lt],(nsuint)1);
            nsuint w=max_threads(ctx.ps[PS_TAIL]);
            if(w>1024) w=1024;
            const mtl_size grid={1,1,1}, group={w,1,1};
            msg<void>(e,RMG_SEL("dispatchThreadgroups:threadsPerThreadgroup:"),grid,group);
            proftail=lt;
            if(profile) flush_profile(PS_TAIL);
        }

        void k_restrict(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            const rmg_kern::rmg_lp F=lp(l);
            obj e=begin(PS_RESTRICT,l);
            bytes_at(e,&T,sizeof(T),0);
            bytes_at(e,&F,sizeof(F),1);
            const int af[]={RMG_F,RMG_U,RMG_P,RMG_N,RMG_S,RMG_W,RMG_E,RMG_T,RMG_B};
            const int ac[]={RMG_F,RMG_U};
            bind(e,l,af,9,2);
            bind(e,l+1,ac,2,11);
            dispatch3(e,PS_RESTRICT,T.chh,2*T.cnx,T.nz);
        }

        void k_prolong(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            obj e=begin(PS_PROLONG,l);
            bytes_at(e,&T,sizeof(T),0);
            const int ac[]={RMG_U};
            const int af[]={RMG_U,RMG_ACT};
            bind(e,l+1,ac,1,2);
            bind(e,l,af,2,3);
            dispatch3(e,PS_PROLONG,T.fhh,2*T.fnx,T.nz);
        }

        void k_coarsen(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            obj e=begin(PS_COARSEN,l);
            bytes_at(e,&T,sizeof(T),0);
            const int af[]={RMG_P,RMG_N,RMG_S,RMG_W,RMG_E,RMG_T,RMG_B,RMG_ACT,RMG_HX,RMG_HY};
            bind(e,l,af,10,2);
            bind(e,l+1,af,10,12);
            dispatch3(e,PS_COARSEN,T.chh,2*T.cnx,T.nz);
        }

        void k_factor(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            if(h16)
            {
                obj e=begin(PS_FACTOR_H,l);
                bytes_at(e,&P,sizeof(P),0);
                const int a[]={RMG_P,RMG_N,RMG_S,RMG_W,RMG_E,RMG_T,RMG_B,RMG_TI,
                               RMG_HN,RMG_HS,RMG_HW,RMG_HE,RMG_HB,RMG_HC,RMG_ACT};
                bind(e,l,a,15,2);
                dispatch3(e,PS_FACTOR_H,P.hh,2*P.nx,1);
                return;
            }
            obj e=begin(PS_FACTOR,l);
            bytes_at(e,&P,sizeof(P),0);
            const int a[]={RMG_P,RMG_T,RMG_B,RMG_TC,RMG_TI,RMG_ACT};
            bind(e,l,a,6,2);
            dispatch3(e,PS_FACTOR,P.hh,2*P.nx,1);
        }

        void k_wrap(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            obj e=begin(PS_WRAP,l);
            bytes_at(e,&P,sizeof(P),0);
            const int a[]={RMG_U};
            bind(e,l,a,1,2);
            dispatch1(e,PS_WRAP,rmg_kern::rmg_nwrap(P));
        }

        void k_zero(int l,int v) override
        {
            obj e=begin(PS_ZERO,l);
            const int a[]={v};
            bind(e,l,a,1,2);
            dispatch1(e,PS_ZERO,lay[l].sk*desc[l].nz);
        }

    private:

        static void* contents(obj b)
        {
            return msg<void*>(b,RMG_SEL("contents"));
        }

        void free_buffers()
        {
            if(!buf.empty()) release(buf[0]);
            release(blk);
            buf.clear();
            blk=nullptr;
        }

        obj new_buffer(long len,std::string &err)
        {
            obj b=msg<obj>(ctx.device,RMG_SEL("newBufferWithLength:options:"),(nsuint)len,MTL_STORAGE_SHARED);
            if(!b)
            {
                err="Metal buffer allocation of "+std::to_string(len)+" bytes failed";
                return nullptr;
            }
            std::memset(contents(b),0,(size_t)len);
            return b;
        }

        //  the current compute encoder, set to kernel k (of level l, for the
        //  profile)
        obj begin(int k,int l)
        {
            curlev=l;
            if(!enc)
            {
                pool p;
                cb=retain(msg<obj>(ctx.queue,RMG_SEL("commandBuffer")));
                enc=retain(msg<obj>(cb,RMG_SEL("computeCommandEncoder")));
            }
            msg<void>(enc,RMG_SEL("setComputePipelineState:"),ctx.ps[k]);
            return enc;
        }

        static void bytes_at(obj e,const void *p,size_t n,int index)
        {
            msg<void>(e,RMG_SEL("setBytes:length:atIndex:"),p,(nsuint)n,(nsuint)index);
        }

        //  arrays a[0..n-1] of level l at buffer indices first, first+1, ...
        void bind(obj e,int l,const int *a,int n,int first)
        {
            obj b[16];
            nsuint o[16];
            for(int m=0;m<n;++m){b[m]=buf[l]; o[m]=(nsuint)(lbase[l]+lay[l].off[a[m]]);}
            const ns_range r={(nsuint)first,(nsuint)n};
            msg<void>(e,RMG_SEL("setBuffers:offsets:withRange:"),(const obj*)b,(const nsuint*)o,r);
        }

        void dispatch1(obj e,int k,long n)
        {
            if(n<=0) return;
            nsuint w=max_threads(ctx.ps[k]);
            if(w>256) w=256;
            const mtl_size grid={(nsuint)n,1,1}, group={w,1,1};
            msg<void>(e,RMG_SEL("dispatchThreads:threadsPerThreadgroup:"),grid,group);
            if(profile) flush_profile(k);
        }

        //  a grid (nx_, ny_, nz_) of work items: groups of up to 64 along the
        //  first, contiguous, dimension and 256 work items in all
        void dispatch3(obj e,int k,long nx_,long ny_,long nz_)
        {
            if(nx_<=0 || ny_<=0 || nz_<=0) return;
            nsuint w=max_threads(ctx.ps[k]);
            if(w>256) w=256;
            nsuint tx=(nsuint)(nx_<64? nx_ : 64);
            nsuint ty=w/tx; if(ty<1) ty=1; if(ty>(nsuint)ny_) ty=(nsuint)ny_;
            nsuint tz=w/(tx*ty); if(tz<1) tz=1; if(tz>(nsuint)nz_) tz=(nsuint)nz_;
            const mtl_size grid={(nsuint)nx_,(nsuint)ny_,(nsuint)nz_}, group={tx,ty,tz};
            msg<void>(e,RMG_SEL("dispatchThreads:threadsPerThreadgroup:"),grid,group);
            if(profile) flush_profile(k);
        }

        //  end the encoder, commit, wait; the GPU time and an error of the buffer
        double commit_and_wait()
        {
            pool p;
            msg<void>(enc,RMG_SEL("endEncoding"));
            msg<void>(cb,RMG_SEL("commit"));
            msg<void>(cb,RMG_SEL("waitUntilCompleted"));
            const double t=msg<double>(cb,RMG_SEL("GPUEndTime"))-msg<double>(cb,RMG_SEL("GPUStartTime"));
            gpu+=t;

            if(msg<long>(cb,RMG_SEL("status"))==MTL_CB_STATUS_ERROR)
            {
                fail="command buffer failed: "+describe(msg<obj>(cb,RMG_SEL("error")));
                failed=true;
            }
            release(enc);
            release(cb);
            enc=nullptr;
            cb=nullptr;
            return t;
        }

        //  profiling: this kernel alone in a command buffer, timed on the GPU
        void flush_profile(int k)
        {
            const double t=commit_and_wait();
            prof[k]+=t;
            ++profn[k];
            if(curlev>=0 && curlev<(int)profl.size()) profl[curlev]+=t;
        }

        metal_ctx &ctx;
        std::vector<obj> buf;      // buffer of every level: level 0's, then blk
        obj blk=nullptr;           // the coarse levels
        obj cb=nullptr;
        obj enc=nullptr;
        long bytes=0;
        bool failed=false;
        std::string fail;
        double gpu=0.0;
        bool profile=false;
        double prof[PS_N];
        long profn[PS_N];
        std::vector<double> profl;
        int curlev=-1, proftail=-1;
    };

    #undef RMG_SEL
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

#endif
