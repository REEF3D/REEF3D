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

#ifndef FIELD_BASE_H_
#define FIELD_BASE_H_

#include "lexer.h"
#include "ckpt_registry.h"

#include <cstddef>
#include <type_traits>

template<typename T>
class field_base
{
public:
    explicit field_base(lexer *p, bool allocate=true) : field_base(p, p->kmax, 0, allocate) {}

    field_base(const field_base&) = delete;
    field_base &operator=(const field_base&) = delete;
    field_base(field_base&&) = delete;
    field_base &operator=(field_base&&) = delete;

    virtual ~field_base()
    {
        ckpt_registry::remove(ckid);

        delete [] V;
        V = nullptr;
    }

    inline T& operator()(int ii, int jj, int kk) noexcept {return V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin];};
    inline T& operator[](int n) noexcept {return V[n];};

    T *data() noexcept {return V;}
    const T *data() const noexcept {return V;}

protected:
    // Vertical-extent-parameterised constructor. operator() folds kz into both
    // the j- and k-strides, so a field with a different number of z-planes needs
    // nothing but a different kz: passing p->kmaxF reproduces the FIJK addressing
    // in iterators3D.h exactly. slack is extra trailing elements, for layouts
    // whose forward-stencil macros reach past the last in-stride slot. See field7.
    field_base(lexer* p, int kz, std::size_t slack, bool allocate=true) :
        imin(p->imin), jkmax(p->jmax*kz), jmin(p->jmin), kmin(p->kmin), kmax(kz)
    {
        size_t _size = static_cast<std::size_t>(p->imax)*jkmax + slack;
        V = allocate ? new T[_size] {} : nullptr;

        if(V!=nullptr)
        {
            if(kmax==p->kmax && slack==0)
            {
                ckid = ckpt_registry::add(p,'C',std::is_same<T,int>::value,V);
            }
            else if(kmax==p->kmaxF)
            {
                // Placeholder for future use, e.g. for a field7 with kmaxF
                // ckid = ckpt_registry::add(p,'F',std::is_same<T,int>::value,V);
            }
        }
    }

    T *V;

private:
    const int imin,jkmax,jmin,kmin,kmax;

    unsigned long long ckid = 0;    // checkpoints (ckpt_registry.h)
};

#endif
