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
        lay.push_back(rmg_layout(lv[l].nx,lv[l].ny,lv[l].nz,h16));
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
        if(per && !z) k_wrap(l);
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
    if(perx || pery) k_wrap(l);
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
        float* F(int l,int a){return (float*)arr(l,a);}
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

        void k_wrap(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int n=rmg_kern::rmg_nwrap(P);
            for(int t=0;t<n;++t)
            rmg_kern::rmg_wrap(P,(unsigned)t,F(l,RMG_U));
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
