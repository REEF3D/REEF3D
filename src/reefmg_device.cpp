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

#include <cstdlib>
#include <cstring>

//  ----------------------------------------------- V-cycle, every device

void reefmg_device::coarsen(int nuse)
{
    for(int l=0;l+1<nuse;++l)
    k_coarsen(l);

    for(int l=0;l<nuse;++l)
    k_factor(l);
}

//  line_zebra_t: colour 0 then 1, or 1 then 0 for dir 1; the periodic halo
//  before every half-sweep (the exchange of reefmg_core::halo)
void reefmg_device::smooth(int l,int sweeps,int dir)
{
    const bool per=(perx || pery);
    for(int sw=0;sw<sweeps;++sw)
    for(int h=0;h<2;++h)
    {
        if(per) k_wrap(l);
        k_smooth(l,(dir==1)? 1-h : h);
    }
}

void reefmg_device::k_coarse_solve(int l,int sweeps)
{
    smooth(l,sweeps,2);
}

void reefmg_device::cycle(int l,int nuse,int pre,int post,int sweepstyle,int coarse_sweeps)
{
    if(l==nuse-1)
    {
        k_coarse_solve(l,coarse_sweeps);
        return;
    }

    const int dpre =(sweepstyle==1)?2:0;
    const int dpost=(sweepstyle==1)?2:1;

    smooth(l,pre,dpre);
    if(perx || pery) k_wrap(l);
    k_residual(l);
    k_restrict(l);
    cycle(l+1,nuse,pre,post,sweepstyle,coarse_sweeps);
    k_prolong(l);
    smooth(l,post,dpost);
}

void reefmg_device::vcycle(int nuse,int pre,int post,int sweepstyle,int coarse_sweeps)
{
    k_zero(0,RMG_U);
    cycle(0,nuse,pre,post,sweepstyle,coarse_sweeps);
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
            for(size_t l=0;l<own.size();++l) std::free(own[l]);
        }

        std::string name() const override {return "reference device (CPU, single threaded)";}

        bool build(const std::vector<rmg_level_desc> &lv,int px,int py,void *ext0,std::string &err) override
        {
            desc=lv; perx=px; pery=py;
            base.assign(lv.size(),0);
            for(size_t l=0;l<lv.size();++l)
            {
                const rmg_layout lay(lv[l].nx,lv[l].ny,lv[l].nz);
                if((long)(lv[l].nx+2)*(lv[l].ny+2)*lv[l].nz>=2147483647L)
                {
                    err="level too large for 32-bit indexing";
                    return false;
                }
                if(l==0 && ext0)
                base[l]=(char*)ext0;
                else
                {
                    base[l]=(char*)std::calloc((size_t)lay.bytes,1);
                    if(!base[l]){err="out of memory"; return false;}
                    own.push_back(base[l]);
                }
            }
            return true;
        }

        void* arr(int l,int a) override
        {
            const rmg_layout lay(desc[l].nx,desc[l].ny,desc[l].nz);
            return base[l]+lay.off[a];
        }

        bool sync(std::string &) override {return true;}

        long memory_bytes() const override
        {
            long s=0;
            for(size_t l=0;l<desc.size();++l)
            if(!(l==0 && own.size()<desc.size()))
            s+=rmg_layout(desc[l].nx,desc[l].ny,desc[l].nz).bytes;
            return s;
        }

    protected:
        float* F(int l,int a){return (float*)arr(l,a);}
        unsigned char* A(int l){return (unsigned char*)arr(l,RMG_ACT);}

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

        void k_smooth(int l,int colour) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int n=rmg_kern::rmg_ncolour(P);
            for(int t=0;t<n;++t)
            rmg_kern::rmg_smooth(P,colour,(unsigned)t,F(l,RMG_F),F(l,RMG_U),F(l,RMG_N),F(l,RMG_S),
                                 F(l,RMG_W),F(l,RMG_E),F(l,RMG_B),F(l,RMG_TC),F(l,RMG_TI),A(l));
        }

        void k_residual(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int n=P.nx*P.ny*P.nz;
            for(int t=0;t<n;++t)
            rmg_kern::rmg_residual(P,(unsigned)t,F(l,RMG_F),F(l,RMG_U),F(l,RMG_R),F(l,RMG_P),F(l,RMG_N),
                                   F(l,RMG_S),F(l,RMG_W),F(l,RMG_E),F(l,RMG_T),F(l,RMG_B));
        }

        void k_restrict(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            const int n=T.cnx*T.cny*T.nz;
            for(int t=0;t<n;++t)
            rmg_kern::rmg_restrict(T,(unsigned)t,F(l,RMG_R),F(l+1,RMG_F),F(l+1,RMG_U));
        }

        void k_prolong(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            const int n=T.fnx*T.fny*T.nz;
            for(int t=0;t<n;++t)
            rmg_kern::rmg_prolong(T,(unsigned)t,F(l+1,RMG_U),F(l,RMG_U),A(l));
        }

        void k_coarsen(int l) override
        {
            const rmg_kern::rmg_tp T=tp(l);
            const int n=T.cnx*T.cny*T.nz;
            for(int t=0;t<n;++t)
            rmg_kern::rmg_coarsen(T,(unsigned)t,
                                  F(l,RMG_P),F(l,RMG_N),F(l,RMG_S),F(l,RMG_W),F(l,RMG_E),F(l,RMG_T),F(l,RMG_B),
                                  A(l),F(l,RMG_HX),F(l,RMG_HY),
                                  F(l+1,RMG_P),F(l+1,RMG_N),F(l+1,RMG_S),F(l+1,RMG_W),F(l+1,RMG_E),
                                  F(l+1,RMG_T),F(l+1,RMG_B),A(l+1),F(l+1,RMG_HX),F(l+1,RMG_HY));
        }

        void k_factor(int l) override
        {
            const rmg_kern::rmg_lp P=lp(l);
            const int n=P.nx*P.ny;
            for(int t=0;t<n;++t)
            rmg_kern::rmg_factor(P,(unsigned)t,F(l,RMG_P),F(l,RMG_T),F(l,RMG_B),F(l,RMG_TC),F(l,RMG_TI),A(l));
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
            const rmg_level_desc &d=desc[l];
            const int n=(d.nx+2)*(d.ny+2)*d.nz;
            float *p=F(l,a);
            for(int t=0;t<n;++t)
            rmg_kern::rmg_zero((unsigned)t,p);
        }

    private:
        std::vector<char*> base;
        std::vector<char*> own;
    };
}

reefmg_device* reefmg_device_create(int kind,std::string &err)
{
    if(kind==2)
    return new reefmg_device_ref();

    if(kind==1)
    {
#ifdef REEF3D_USE_METAL
        return reefmg_device_metal_create(err);
#else
        err="this build has no Metal GPU support (build with make METAL=1 or cmake -DREEF3D_USE_METAL=ON on macOS)";
        return 0;
#endif
    }

    err="unknown device "+std::to_string(kind);
    return 0;
}
