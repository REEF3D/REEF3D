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
    {cache_addressing();}

    slice_base(const slice_base&) = delete;
    slice_base& operator=(const slice_base&) = delete;
    slice_base(slice_base&&) = delete;
    slice_base& operator=(slice_base&&) = delete;

    virtual ~slice_base() = default;

    /// Width of the cached stride. Long for every payload — see
    /// field_base::stride_t for the measured rationale (chiefly: every
    /// container a loop body touches must agree on the width, or their index
    /// computations cannot be shared).
    using stride_t = long;

    inline T& operator()(int ii, int jj) noexcept {return m_base[ii*m_js + jj];}
    inline const T& operator()(int ii, int jj) const noexcept {return m_base[ii*m_js + jj];}

    /*!
     * @brief Lightweight, capture-by-value view: the folded base pointer and
     * stride copied out once, so a loop body indexes local values instead of
     * reloading this->m_base/m_js through operator(). Mirrors the amrex::Array4
     * captured by FIELDLOOP under USE_AMREX (materialized once per tile, then
     * indexed by value) so the same LOOP macros can be hoisted out of the loop
     * body and stay vectorizable.
     */
    struct View
    {
        T *__restrict base;        ///< already origin-folded, see cache_addressing()
        stride_t js;
        inline T &operator()(int ii, int jj) const noexcept
        { return base[ii*js + jj]; }
    };

    struct ConstView
    {
        const T *__restrict base;
        stride_t js;
        inline const T &operator()(int ii, int jj) const noexcept
        { return base[ii*js + jj]; }
    };

    inline View view() noexcept { return View{m_base, m_js}; }
    inline ConstView const_view() const noexcept { return ConstView{m_base, m_js}; }

    T *data() noexcept {return V.data();}
    const T *data() const noexcept {return V.data();}

    // whole array including ghost cells
    std::size_t size() const noexcept {return V.size();}
    T *begin() noexcept {return V.data();}
    const T *begin() const noexcept {return V.data();}
    T *end() noexcept {return V.data()+V.size();}
    const T *end() const noexcept {return V.data()+V.size();}

protected:
    /*!
     * @brief Precomputes the folded addressing used by the derived operator().
     *
     * 2D counterpart of field_base::cache_addressing — see there for the
     * rationale. Here
     *
     *   V[(ii-imin)*jmax + (jj-jmin)] == m_base[ii*m_js + jj]
     *
     * Must be called after anything that reseats V.
     */
    void cache_addressing() noexcept
    {
        m_js   = jmax;
        m_base = V.data() - static_cast<long>(imin)*jmax - jmin;
    }

    std::vector<T> V;

private:
    T       *m_base = nullptr; ///< origin-folded base: m_base[ii*m_js + jj]
    stride_t m_js   = 0;       ///< i-stride (jmax); long for every payload, see stride_t

    const int imin,jmin;
    const stride_t jmax; ///< i-stride; stride_t like m_js, which it seeds
};

#endif
