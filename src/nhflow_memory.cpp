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


#include"nhflow_memory.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ckpt_registry.h"
#include<sys/mman.h>
#include<unistd.h>
#include<cstdint>
#include<cstring>
#include<vector>
#include<iostream>

namespace
{
// The whole pages of an array that hold zeros only and are in memory are given back to the system (madvise
// MADV_DONTNEED). Such a page of private anonymous memory reads as zeros again afterwards, so nothing changes for the
// run; it takes memory again only when it is written. Returns the bytes given back.
size_t release_zero_pages(const void *ptr, size_t bytes)
{
    if(ptr==nullptr || bytes==0)
    return 0;

    const uintptr_t ps = (uintptr_t)sysconf(_SC_PAGESIZE);
    const uintptr_t a0 = ((uintptr_t)ptr + ps - 1) & ~(ps - 1);
    const uintptr_t e = ((uintptr_t)ptr + bytes) & ~(ps - 1);

    if(e<=a0)
    return 0;

    // only the pages in memory: the others are given back already (or were never touched)
    const size_t np = (e - a0)/ps;
    // mincore takes char* on macOS and unsigned char* on Linux
#if defined(__APPLE__)
    std::vector<char> res(np);
#else
    std::vector<unsigned char> res(np);
#endif
    if(mincore((void*)a0, e - a0, res.data())!=0)
    return 0;

    size_t freed = 0;
    size_t q = 0;

    while(q<np)
    {
        // a run of resident zero pages
        size_t r = q;

        while(r<np && (res[r]&1))
        {
            const uint64_t *w = (const uint64_t*)(a0 + r*ps);
            bool zero = true;

            for(size_t n=0; n<ps/sizeof(uint64_t); ++n)
            if(w[n]!=0)
            {
                zero = false;
                break;
            }

            if(!zero)
            break;

            ++r;
        }

        if(r>q)
        {
            if(madvise((void*)(a0 + q*ps), (r - q)*ps, MADV_DONTNEED)==0)
            freed += (r - q)*ps;

            q = r;
        }
        else
        ++q;
    }

    return freed;
}
}

void nhflow_memory_trim(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    size_t freed = 0;

    // the 3D fields of fdm_nhf
    const size_t n3 = size_t(p->imax)*p->jmax*(p->kmax+2);

    double *f3[] = {d->U,d->V,d->W,d->omegaF,d->UH,d->VH,d->WH,d->P,d->RO,d->VISC,d->EV,d->EV0,d->EVH,
                    d->F,d->G,d->H,d->L,d->Fext,d->Gext,d->Hext,d->POR,d->PORPART,d->PORDEM,d->test,d->KIN,d->CONC,
                    d->SALT,d->TEMP,d->SOLID,d->FB,d->FHB,d->PORSTRUC,d->Fx,d->Fy,d->Fz,d->FEx,d->FEy,
                    d->Ss,d->Sn,d->Se,d->Sw,d->Un,d->Us,d->Ub,d->Ut,d->Ve,d->Vw,d->Vb,d->Vt,d->Wb,d->Wt,
                    d->Vn,d->Vs,d->Wn,d->Ws,d->Ue,d->Uw,d->We,d->Ww,
                    d->UHn,d->UHs,d->UHe,d->UHw,d->VHn,d->VHs,d->VHe,d->VHw,d->WHn,d->WHs,d->WHe,d->WHw};

    for(double *f : f3)
    freed += release_zero_pages(f, n3*sizeof(double));

    // every slice of this grid (the slices register in ckpt_registry: fdm_nhf, the momentum classes, ...)
    const size_t n2 = size_t(p->imax)*p->jmax;

    ckpt_registry::visit(p, [&](const ckpt_registry::entry &e)
    {
        if(e.kind=='S')
        freed += release_zero_pages(e.V, n2*(e.isint ? sizeof(int) : sizeof(double)));
    });

    // the matrix and right-hand side of the implicit 3D solves (diffusion, turbulence)
    for(std::vector<double> *v : {&d->M.n,&d->M.s,&d->M.e,&d->M.w,&d->M.b,&d->M.t,&d->M.p,&d->rhsvec.V})
    freed += release_zero_pages(v->data(), v->size()*sizeof(double));

    const double mb = pgc->globalsum(double(freed))/1.0e6;

    if(p->mpirank==0)
    std::cout<<"NHFLOW semi-implicit: "<<int(mb+0.5)<<" MB of zero pages given back after the initialisation"<<std::endl;
}
