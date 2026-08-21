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

#include <cstddef>
#include <vector>

template<typename T>
class field_base
{
public:
    field_base(lexer *p, bool allocate=true) : field_base(p, p->kmax, 0, allocate) {}

    field_base(const field_base&) = delete;
    field_base& operator=(const field_base&) = delete;
    field_base(field_base&&) = delete;
    field_base& operator=(field_base&&) = delete;

    virtual ~field_base() = default;

    inline T& operator()(int ii, int jj, int kk) noexcept {return V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin];};

    T *data() noexcept {return V.data();}
    const T *data() const noexcept {return V.data();}

    // whole array including ghost cells and slack; 0 for allocate=false
    std::size_t size() const noexcept {return V.size();}
    T *begin() noexcept {return V.data();}
    const T *begin() const noexcept {return V.data();}
    T *end() noexcept {return V.data()+V.size();}
    const T *end() const noexcept {return V.data()+V.size();}

protected:
    // Vertical-extent-parameterised constructor. operator() folds kz into both
    // the j- and k-strides, so a field with a different number of z-planes needs
    // nothing but a different kz: passing p->kmaxF reproduces the FIJK addressing
    // in iterators3D.h exactly. slack is extra trailing elements, for layouts
    // whose forward-stencil macros reach past the last in-stride slot. See field7.
    // allocate=false: the field is not used in this run (V stays empty)
    field_base(lexer *p, int kz, std::size_t slack, bool allocate=true) :
        V(allocate ? static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax)*static_cast<std::size_t>(kz) + slack : 0, T{}),
        imin(p->imin), jmin(p->jmin),
        kmin(p->kmin), kmax(kz),
        jkmax(p->jmax*kz)
    {}

    std::vector<T> V;

private:
    const int imin,jmin,kmin,kmax,jkmax;
};

#endif
