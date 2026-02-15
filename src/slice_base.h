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

#include <cstddef>
#include <type_traits>

template<typename T>
class slice_base
{
public:
    explicit slice_base(lexer *p) : imin(p->imin), jmin(p->jmin), jmax(p->jmax),
        m_size(static_cast<std::size_t>(p->imax)*jmax)
    {
        V = new T[m_size] {};

        ckid = ckpt_registry::add(p,'S',std::is_same<T,int>::value,V);
    }

    slice_base(const slice_base&) = delete;
    slice_base &operator=(const slice_base&) = delete;
    slice_base(slice_base&&) = delete;
    slice_base &operator=(slice_base&&) = delete;

    virtual ~slice_base()
    {
        ckpt_registry::remove(ckid);

        delete [] V;
        V = nullptr;
    }

    inline T& operator()(int ii, int jj) noexcept {return V[(ii-imin)*jmax + (jj-jmin)];};

    T *data() noexcept {return V;}
    const T *data() const noexcept {return V;}

    // whole array including ghost cells
    std::size_t size() const noexcept {return m_size;}
    T *begin() noexcept {return V;}
    const T *begin() const noexcept {return V;}
    T *end() noexcept {return V+m_size;}
    const T *end() const noexcept {return V+m_size;}

protected:
    T *V;

private:
    const int imin,jmin,jmax;
    const std::size_t m_size;

    unsigned long long ckid = 0;    // checkpoints (ckpt_registry.h)
};

#endif
