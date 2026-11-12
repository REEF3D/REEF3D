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

#include "reefmg_device.h"
#include "reefmg_device_kernels.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>

//  ----------------------------------------------- V-cycle, every device

reefmg_device::reefmg_device()
{
    const char *e=std::getenv("REEFMG_DEVICE_TAIL");
    if(e && e[0]) tailcols=std::atol(e);
}

long reefmg_device::plan_levels(const std::vector<rmg_level_desc> &lv,int px,int py,bool half16)
{
    desc=lv; perx=px; pery=py; h16=half16;
    lay.clear();
    lbase.assign(lv.size(),0);
    long o=0;
    for(size_t l=0;l<lv.size();++l)
    {
        lay.push_back(rmg_layout(lv[l].nx,lv[l].ny,lv[l].nz,h16,l==0 && lv[l].kry!=0));
        if(l>0)
        {
            lbase[l]=o;
            o+=(lay[l].bytes+255)/256*256;
        }
    }
    return o;
}

rmg_kern::rmg_lp reefmg_device::lp(int l) const
{
    const rmg_level_desc &d=desc[l];
    rmg_kern::rmg_lp P;
    P.nx=d.nx; P.ny=d.ny; P.nz=d.nz; P.hh=lay[l].hh; P.si=lay[l].si; P.sk=(int)lay[l].sk;
    P.perx=perx; P.pery=pery; P.h16=h16? 1 : 0;
    return P;
}

rmg_kern::rmg_tp reefmg_device::tp(int l) const
{
    const rmg_level_desc &f=desc[l], &c=desc[l+1];
    rmg_kern::rmg_tp T;
    T.fnx=f.nx; T.fny=f.ny; T.nz=f.nz; T.fhh=lay[l].hh; T.fsi=lay[l].si; T.fsk=(int)lay[l].sk;
    T.cnx=c.nx; T.cny=c.ny; T.chh=lay[l+1].hh; T.csi=lay[l+1].si; T.csk=(int)lay[l+1].sk;
    T.rx=c.rx; T.ry=c.ry; T.perx=perx; T.pery=pery;
    return T;
}

int reefmg_device::tail_level(int nuse) const
{
    const int maxlev=12;    // rmg_tail::L
    int lt=nuse;
    if(tailcols<0)
    return lt;
    for(int l=nuse-1;l>=1 && nuse-l<=maxlev;--l)
    {
        const long cols=(long)desc[l].nx*lay[l].hh;
        if(cols<=tailcols && lbase[nuse-1]+lay[nuse-1].bytes-lbase[l]<2147483647L) lt=l;
        else break;
    }
    if(lt==nuse && nuse>=2 && (long)desc[nuse-1].nx*lay[nuse-1].hh<=4096)
    lt=nuse-1;
    return lt;
}

void reefmg_device::tail_params(int lt,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,rmg_kern::rmg_tail &X) const
{
    std::memset(&X,0,sizeof(X));
    X.nlev=nuse-lt;
    X.pre=pre; X.post=post;
    X.dpre =(sweepstyle==1)?2:0;
    X.dpost=(sweepstyle==1)?2:1;
    X.csweeps=coarse_sweeps;
    for(int l=lt;l<nuse;++l)
    {
        rmg_kern::rmg_tlv &V=X.L[l-lt];
        V.P=lp(l);
        if(l+1<nuse) V.T=tp(l);
        const long b=lbase[l]-lbase[lt];
        const long *o=lay[l].off;
        V.of=(int)(b+o[RMG_F]);   V.ou=(int)(b+o[RMG_U]);   V.op=(int)(b+o[RMG_P]);
        V.on=(int)(b+o[RMG_N]);   V.os=(int)(b+o[RMG_S]);   V.ow=(int)(b+o[RMG_W]);
        V.oe=(int)(b+o[RMG_E]);   V.ot=(int)(b+o[RMG_T]);   V.ob=(int)(b+o[RMG_B]);
        V.otc=(int)(b+o[RMG_TC]); V.oti=(int)(b+o[RMG_TI]); V.oact=(int)(b+o[RMG_ACT]);
        V.ohn=(int)(b+o[RMG_HN]); V.ohs=(int)(b+o[RMG_HS]); V.ohw=(int)(b+o[RMG_HW]);
        V.ohe=(int)(b+o[RMG_HE]); V.ohb=(int)(b+o[RMG_HB]); V.ohc=(int)(b+o[RMG_HC]);
    }
}

void reefmg_device::coarsen(int nuse)
{
    for(int l=0;l+1<nuse;++l)
    k_coarsen(l);

    for(int l=0;l<nuse;++l)
    k_factor(l);
}

//  line_zebra_t: colour 0 then 1, or 1 then 0 for dir 1; the periodic halo
//  before every half-sweep (the exchange of reefmg_core::halo).  from_zero:
//  u is to be taken as zero - the first half-sweep neither reads nor needs
//  the other colour and sets it to zero.
void reefmg_device::smooth(int l,int sweeps,int dir,bool from_zero)
{
    const bool per=(perx || pery);
    for(int sw=0;sw<sweeps;++sw)
    for(int h=0;h<2;++h)
    {
        const bool z=(from_zero && sw==0 && h==0);
        if(per && !z) k_wrap(l,RMG_U);
        k_smooth(l,(dir==1)? 1-h : h,z? 1 : 0);
    }
}

//  the tail step by step: the V-cycle from level lt, recursively
void reefmg_device::k_tail(int lt,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps)
{
    const int t=tail0;
    tail0=-1;
    cycle(lt,nuse,pre,post,sweepstyle,coarse_sweeps);
    tail0=t;
}

//  The residual is computed inside the restriction (k_restrict), so it is
//  never stored; the periodic halo of u first.
void reefmg_device::cycle(int l,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,bool from_zero)
{
    if(l==tail0 && l>0)
    {
        k_tail(l,nuse,pre,post,sweepstyle,coarse_sweeps);
        return;
    }
    if(l==nuse-1)
    {
        if(from_zero) k_zero(l,RMG_U);
        smooth(l,coarse_sweeps,2);
        return;
    }

    const int dpre =(sweepstyle==1)?2:0;
    const int dpost=(sweepstyle==1)?2:1;

    if(from_zero && pre==0) k_zero(l,RMG_U);
    smooth(l,pre,dpre,from_zero && pre>0);
    if(perx || pery) k_wrap(l,RMG_U);
    k_restrict(l);
    cycle(l+1,nuse,pre,post,sweepstyle,coarse_sweeps);
    k_prolong(l);
    smooth(l,post,dpost);
}

//  u(0) from zero: the first pre-smoothing half-sweep starts from u = 0 and
//  replaces zeroing it; further cycles start from the result
void reefmg_device::vcycle(int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,int ncycle)
{
    tail0=tail_level(nuse);
    for(int c=0;c<std::max(ncycle,1);++c)
    cycle(0,nuse,pre,post,sweepstyle,coarse_sweeps,c==0);
    tail0=-1;
    nvc+=std::max(ncycle,1);
}

//  ------------------------------------------- Krylov solve, every device

void reefmg_device::kprecond(int x,int y,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,int ncycle)
{
    fa=x; ua=y;
    vcycle(nuse,pre,post,sweepstyle,coarse_sweeps,ncycle);
    fa=RMG_F; ua=RMG_U;
    //  the operator reads the halo of y
    if(perx || pery) k_wrap(0,y);
}

bool reefmg_device::kdots(double &d1,double &d2)
{
    std::string err;
    d1=d2=0.0;
    if(!sync(err))
    return false;
    const float *p=(const float*)arr(0,RMG_KD);
    for(int t=0;t<RMG_DOTG;++t)
    {
        d1+=p[t];
        d2+=p[RMG_DOTG+t];
    }
    return true;
}

//  BiCGStab of reefmg_core::solve in single precision on the device's fine
//  operator, right-preconditioned by the V-cycle: x in RMG_U, s in RMG_KR.
//  The shadow residual is f itself (r0 = f, x0 = 0, and f is not changed),
//  r is first formed as s = f - alpha v and x first set at its first update,
//  and the dot products of the new s and r come with their updates
//  (k_axpy3_dot).  A failed device stops the iteration; the caller's sync
//  reports it.
int reefmg_device::krylov(int nuse,int pre,int post,int sweepstyle,int coarse_sweeps,int ncycle,
                          double tol,int maxit,double &relres)
{
    const int R=RMG_KR, H=RMG_F, P=RMG_KP, V=RMG_KV, Y=RMG_KY, Z=RMG_KZ, T=RMG_KT;
    relres=1.0;
    bool xset=false;    // x (RMG_U) written

    k_dot2(0,RMG_F,RMG_F,RMG_F,RMG_F);
    double bb,dum;
    if(!kdots(bb,dum)) return 0;
    if(!(bb>0.0))
    {
        k_axpy3(0,RMG_U,0.0f,RMG_F,0.0f,RMG_F,0.0f);
        relres=0.0;
        return 0;
    }
    const double bn=std::sqrt(bb);

    double rho=1.0,alpha=1.0,omega=1.0,rho1=bb,rn=bn;
    int it=0;

    while(rn/bn>tol && it<maxit)
    {
        if(std::fabs(rho1)<1.0e-300) break;

        //  p = r, or r + beta (p - omega v); r = f at the start
        if(it==0)
        k_axpy3(0,P,0.0f,RMG_F,1.0f,RMG_F,0.0f);
        else
        {
            const double beta=(rho1/rho)*(alpha/omega);
            k_axpy3(0,P,(float)beta,R,1.0f,V,(float)(-beta*omega));
        }
        rho=rho1;

        kprecond(P,Y,nuse,pre,post,sweepstyle,coarse_sweeps,ncycle);
        k_apply(0,-1,Y,V);
        k_dot2(0,H,V,H,V);
        double den;
        if(!kdots(den,dum)) break;
        if(std::fabs(den)<1.0e-300) break;
        alpha=rho/den;

        //  s = r - alpha v, |s|^2
        if(it==0)
        k_axpy3_dot(0,R,0.0f,RMG_F,1.0f,V,(float)(-alpha),R);
        else
        k_axpy3_dot(0,R,1.0f,V,(float)(-alpha),V,0.0f,R);
        double ss;
        if(!kdots(ss,dum)) break;
        const double sn=std::sqrt(std::max(ss,0.0));
        if(sn/bn<=tol)
        {
            k_axpy3(0,RMG_U,xset? 1.0f : 0.0f,Y,(float)alpha,Y,0.0f);
            xset=true;
            rn=sn;
            ++it;
            break;
        }

        kprecond(R,Z,nuse,pre,post,sweepstyle,coarse_sweeps,ncycle);
        k_apply(0,-1,Z,T);
        k_dot2(0,T,T,T,R);
        double tt,ts;
        if(!kdots(tt,ts)) break;
        omega=(tt>0.0)? ts/tt : 0.0;

        //  x += alpha y + omega z;  r = s - omega t, |r|^2 and (rhat,r)
        k_axpy3(0,RMG_U,xset? 1.0f : 0.0f,Y,(float)alpha,Z,(float)omega);
        xset=true;
        k_axpy3_dot(0,R,1.0f,T,(float)(-omega),T,0.0f,H);
        double rr;
        if(!kdots(rr,rho1)) break;
        rn=std::sqrt(std::max(rr,0.0));
        ++it;

        if(std::fabs(omega)<1.0e-300) break;
    }

    if(!xset)
    k_axpy3(0,RMG_U,0.0f,RMG_F,0.0f,RMG_F,0.0f);

    relres=rn/bn;
    return it;
}

//  ------------------------------------------------------ reference device
//  The kernels of reefmg_device_kernels.h in plain loops, one work item
//  after the other: exactly the code the GPU executes, in single precision,
//  so it tests the device path - layout, transfers, coarsening, the host side
//  in reefmg_core - on any machine.  Not meant to be fast.
namespace
{
    class reefmg_device_ref final : public reefmg_device
    {
    public:
        ~reefmg_device_ref() override
        {
            release();
        }

        std::string name() const override {return "reference device (CPU, single threaded)";}

        bool build(const std::vector<rmg_level_desc> &lv,int px,int py,bool half16,void *ext0,std::string &err) override
        {
            release();
            const long nb=plan_levels(lv,px,py,half16);
            for(size_t l=0;l<lv.size();++l)
            if(lay[l].sk*lv[l].nz>=2147483647L)
            {
                err="level too large for 32-bit indexing";
                return false;
            }
            base.assign(lv.size(),0);
            if(ext0)
            base[0]=(char*)ext0;
            else
            {
                own0=(char*)std::calloc((size_t)lay[0].bytes,1);
                if(!own0){err="out of memory"; return false;}
                base[0]=own0;
            }
            if(nb>0)
            {
                block=(char*)std::calloc((size_t)nb,1);
                if(!block){err="out of memory"; return false;}
                for(size_t l=1;l<lv.size();++l) base[l]=block+lbase[l];
            }
            blockbytes=nb;
            return true;
        }

        void* arr(int l,int a) override
        {
            return base[l]+lay[l].off[a];
        }

        bool sync(std::string &) override {return true;}

        long memory_bytes() const override
        {
            return blockbytes+(own0? lay[0].bytes : 0);
        }

    protected:
        //  the arrays of the kernels: f and u of level 0 through amap()
        float* F(int l,int a){return (float*)arr(l,amap(l,a));}
        rmg_kern::rmg_h* H(int l,int a){return (rmg_kern::rmg_h*)arr(l,a);}
        unsigned char* A(int l){return (unsigned char*)arr(l,RMG_ACT);}

        //  the grids of reefmg_device_kernels.h: (m, row-colour, layer)
        void k_smooth(int l,int colour,int zero) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            for(int i=0;i<P.nx;++i)
            for(int m=0;m<P.hh;++m)
            if(h16)
            rmg_kern::rmg_smooth_h(P,colour,m,i,zero,F(l,RMG_F),F(l,RMG_U),H(l,RMG_HN),H(l,RMG_HS),
                                   H(l,RMG_HW),H(l,RMG_HE),H(l,RMG_HB),H(l,RMG_HC),F(l,RMG_TI),A(l));
            else
            rmg_kern::rmg_smooth(P,colour,m,i,zero,F(l,RMG_F),F(l,RMG_U),F(l,RMG_N),F(l,RMG_S),
                                 F(l,RMG_W),F(l,RMG_E),F(l,RMG_B),F(l,RMG_TC),F(l,RMG_TI),A(l));
        }

        //  the tail kernel's code, one work item doing all of it
        void k_tail(int lt,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps) override
        {
            rmg_kern::rmg_tail X;
            tail_params(lt,nuse,pre,post,sweepstyle,coarse_sweeps,X);
            rmg_kern::rmg_tail_run(X,base[lt],0u,1u);
        }

        void k_restrict(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            const rmg_kern::rmg_lp P=lp(l);
            for(int k=0;k<T.nz;++k)
            for(int r=0;r<2*T.cnx;++r)
            for(int m=0;m<T.chh;++m)
            rmg_kern::rmg_restrict_residual(T,P,m,r,k,F(l,RMG_F),F(l,RMG_U),F(l,RMG_P),F(l,RMG_N),F(l,RMG_S),
                                            F(l,RMG_W),F(l,RMG_E),F(l,RMG_T),F(l,RMG_B),F(l+1,RMG_F),F(l+1,RMG_U));
        }

        void k_prolong(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            for(int k=0;k<T.nz;++k)
            for(int r=0;r<2*T.fnx;++r)
            for(int m=0;m<T.fhh;++m)
            rmg_kern::rmg_prolong(T,m,r,k,F(l+1,RMG_U),F(l,RMG_U),A(l));
        }

        void k_coarsen(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            for(int k=0;k<T.nz;++k)
            for(int r=0;r<2*T.cnx;++r)
            for(int m=0;m<T.chh;++m)
            rmg_kern::rmg_coarsen(T,m,r,k,
                                  F(l,RMG_P),F(l,RMG_N),F(l,RMG_S),F(l,RMG_W),F(l,RMG_E),F(l,RMG_T),F(l,RMG_B),
                                  A(l),F(l,RMG_HX),F(l,RMG_HY),
                                  F(l+1,RMG_P),F(l+1,RMG_N),F(l+1,RMG_S),F(l+1,RMG_W),F(l+1,RMG_E),
                                  F(l+1,RMG_T),F(l+1,RMG_B),A(l+1),F(l+1,RMG_HX),F(l+1,RMG_HY));
        }

        void k_factor(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            for(int r=0;r<2*P.nx;++r)
            for(int m=0;m<P.hh;++m)
            if(h16)
            rmg_kern::rmg_factor_h(P,m,r,F(l,RMG_P),F(l,RMG_N),F(l,RMG_S),F(l,RMG_W),F(l,RMG_E),
                                   F(l,RMG_T),F(l,RMG_B),F(l,RMG_TI),H(l,RMG_HN),H(l,RMG_HS),
                                   H(l,RMG_HW),H(l,RMG_HE),H(l,RMG_HB),H(l,RMG_HC),A(l));
            else
            rmg_kern::rmg_factor(P,m,r,F(l,RMG_P),F(l,RMG_T),F(l,RMG_B),F(l,RMG_TC),F(l,RMG_TI),A(l));
        }

        void k_wrap(int l,int a) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int n=rmg_kern::rmg_nwrap(P);
            for(int t=0;t<n;++t)
            rmg_kern::rmg_wrap(P,(unsigned)t,F(l,a));
        }

        void k_apply(int l,int f,int u,int y) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int mode=(f>=0)? 1 : 0;
            float *ff=F(l,(f>=0)? f : u);
            for(int k=0;k<P.nz;++k)
            for(int r=0;r<2*P.nx;++r)
            for(int m=0;m<P.hh;++m)
            rmg_kern::rmg_apply(P,m,r,k,mode,ff,F(l,u),F(l,y),F(l,RMG_P),F(l,RMG_N),F(l,RMG_S),
                                F(l,RMG_W),F(l,RMG_E),F(l,RMG_T),F(l,RMG_B));
        }

        void k_axpy3(int l,int w,float c0,int a,float c1,int b,float c2) override
        {
            rmg_kern::rmg_vp V;
            V.c0=c0; V.c1=c1; V.c2=c2; V.n=(int)(lay[l].sk*desc[l].nz); V.g=RMG_DOTG;
            float *pw=F(l,w);
            const float *pa=F(l,a), *pb=F(l,b);
            for(int t=0;t<V.n;++t)
            rmg_kern::rmg_axpy3(V,(unsigned)t,pw,pa,pb);
        }

        void k_dot2(int l,int a,int b,int c,int d) override
        {
            rmg_kern::rmg_vp V;
            V.c0=V.c1=V.c2=0.0f; V.n=(int)(lay[l].sk*desc[l].nz); V.g=RMG_DOTG;
            float *part=F(l,RMG_KD);
            for(int t=0;t<V.g;++t)
            rmg_kern::rmg_dot2(V,(unsigned)t,F(l,a),F(l,b),F(l,c),F(l,d),part);
        }

        void k_axpy3_dot(int l,int w,float c0,int a,float c1,int b,float c2,int h) override
        {
            rmg_kern::rmg_vp V;
            V.c0=c0; V.c1=c1; V.c2=c2; V.n=(int)(lay[l].sk*desc[l].nz); V.g=RMG_DOTG;
            float *part=F(l,RMG_KD), *pw=F(l,w);
            for(int t=0;t<V.g;++t)
            rmg_kern::rmg_axpy3_dot(V,(unsigned)t,pw,F(l,a),F(l,b),F(l,h),part);
        }

        void k_zero(int l,int a) override
        {
            const int n=(int)(lay[l].sk*desc[l].nz);
            float *p=F(l,a);
            for(int t=0;t<n;++t)
            rmg_kern::rmg_zero((unsigned)t,p);
        }

    private:
        void release()
        {
            std::free(own0); own0=0;
            std::free(block); block=0;
            blockbytes=0;
        }

        std::vector<char*> base;
        char *own0=0, *block=0;
        long blockbytes=0;
    };
}

reefmg_device* reefmg_device_create(int kind,std::string &err)
{
    if(kind==2)
    return new reefmg_device_ref();

    if(kind==1)
    {
#if REEF3D_METAL
        return reefmg_device_metal_create(err);
#else
        err="no Metal GPU support in this build (macOS only)";
        return 0;
#endif
    }

    err="unknown device "+std::to_string(kind);
    return 0;
}
