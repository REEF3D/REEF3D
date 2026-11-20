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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef SLICE_BASE_H_
#define SLICE_BASE_H_

#include "lexer.h"
#include "ckpt_registry.h"
#include <type_traits>
#include <cstdlib>
#include <new>
#ifdef NHF_KCHECK
#include"nhf_kcheck.h"
#endif

template<typename T>
class slice_base
{
public:
    // the zeros of a new slice from calloc: its pages take memory only when written (a slice an option never uses
    // costs nothing)
    slice_base(lexer* p) : imin(p->imin), imax(p->imax), jmin(p->jmin), jmax(p->jmax)
    {
        V = zeros();
        ckid = ckpt_registry::add(p,'S',std::is_same<T,int>::value,V);
    }
    // on = false: the slice of an option that is not used; no array (V = nullptr), not in the checkpoints
    slice_base(lexer* p, bool on) : imin(p->imin), imax(p->imax), jmin(p->jmin), jmax(p->jmax)
    {
        V = on ? zeros() : nullptr;
        if(on)
        ckid = ckpt_registry::add(p,'S',std::is_same<T,int>::value,V);
    }
    virtual ~slice_base()
    {
        ckpt_registry::remove(ckid);
        free(V);
        V = nullptr;
    }

    slice_base(const slice_base&) = delete;
    slice_base& operator=(const slice_base&) = delete;
    slice_base(slice_base&&) = delete;
    slice_base& operator=(slice_base&&) = delete;

#ifdef NHF_KCHECK
    // check build (nhf_kcheck.h): a row outside the slice is reported with the address of the caller
    __attribute__((noinline)) T& operator()(int ii, int jj) noexcept
    {
        if(jmin>-3 && (jj<jmin || jj>=jmin+jmax))
        nhf_jslice(__builtin_return_address(0),jj,jmin,jmax);
        return V[(ii-imin)*jmax + (jj-jmin)];
    }
#else
    inline T& operator()(int ii, int jj) noexcept {return V[(ii-imin)*jmax + (jj-jmin)];};
#endif
    inline T& operator[](int n) noexcept {return V[n];};

    T *V;
protected:
    T* zeros() const
    {
        T *v = static_cast<T*>(calloc(size_t(imax)*size_t(jmax),sizeof(T)));
        if(v==nullptr)
        throw std::bad_alloc();
        return v;
    }

    const int imin,imax,jmin,jmax;
    unsigned long long ckid = 0;    // checkpoints (ckpt_registry.h)
};

#endif
