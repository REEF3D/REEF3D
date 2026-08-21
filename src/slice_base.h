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

#include <cstddef>
#include <vector>

template<typename T>
class slice_base
{
public:
    slice_base(lexer *p) :
        V(static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax), T{}),
        imin(p->imin), jmin(p->jmin),
        jmax(p->jmax)
    {};

    slice_base(const slice_base&) = delete;
    slice_base& operator=(const slice_base&) = delete;
    slice_base(slice_base&&) = delete;
    slice_base& operator=(slice_base&&) = delete;

    virtual ~slice_base() = default;

    inline T& operator()(int ii, int jj) noexcept {return V[(ii-imin)*jmax + (jj-jmin)];};

    T *data() noexcept {return V.data();}
    const T *data() const noexcept {return V.data();}

    // whole array including ghost cells
    std::size_t size() const noexcept {return V.size();}
    T *begin() noexcept {return V.data();}
    const T *begin() const noexcept {return V.data();}
    T *end() noexcept {return V.data()+V.size();}
    const T *end() const noexcept {return V.data()+V.size();}

protected:
    std::vector<T> V;

private:
    const int imin,jmin,jmax;
    const std::size_t n;
};

#endif
