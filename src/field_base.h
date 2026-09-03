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
#include "looping.h"

#include <algorithm>
#include <cstddef>
#include <vector>

template<typename T>
class field_base
{
public:
    field_base(lexer *pp, bool allocate=true) : field_base(pp, pp->kmax, 0, allocate) {}

    field_base(const field_base&) = delete;
    field_base& operator=(const field_base&) = delete;
    field_base(field_base&&) = delete;
    field_base& operator=(field_base&&) = delete;

    virtual ~field_base() = default;

    /*!
     * @brief Width of the cached strides. Long for every payload, deliberately.
     *
     * Two reasons, both measured:
     *
     *  1. An int store through the reference operator() returns may alias an
     *     int stride member under TBAA, forcing the strides to be reloaded on
     *     every iteration of a writing loop. Long removes that alias and lets
     *     the address strength-reduce to an induction variable (~1.7-2.8x on
     *     int payloads).
     *
     *  2. More importantly, the width must be the SAME across every container a
     *     loop body touches. A REEF3D loop is typically "test an int flag,
     *     update a double field"; if the two disagree on stride width their
     *     index computations cannot be shared and the body builds two of them.
     *     Uniform long measured 1.35-1.50x on that shape versus 0.99-1.10x for
     *     uniform int and 1.15-1.27x for a mixed int/long split.
     *
     * The one shape narrow strides win (a pure double constant-store, ~1.2x) is
     * essentially just setVal, which takes the std::fill path anyway.
     * Do not make this payload-dependent — that was tried and measured worse.
     */
    using stride_t = long;

    // Origin and both strides are folded into m_base by cache_addressing(), so
    // these reach the element with two multiplies and no member subtractions.
    // Equivalent to V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin].
    inline T& operator()(int ii, int jj, int kk) noexcept {return m_base[ii*m_js + jj*m_ks + kk];};
    inline const T& operator()(int ii, int jj, int kk) const noexcept {return m_base[ii*m_js + jj*m_ks + kk];};

    // includeGhost=false: LOOP cells only (PCHECK: flag4 > 0). Cells with
    // flag4 <= 0 (solid, object, in/outflow, ...) keep their values, as do the
    // ghosts. Same as ArrayWrapper3D.
    void setVal(T val, bool includeGhost = false)
    {
        if(V.empty())
        {
            return; // allocate=false: no storage to fill
        }
        else if(includeGhost)
        {
            std::fill(V.begin(), V.end(), val);
        }
        else
        {
            int i,j,k;
            LOOP
            {
                operator()(i,j,k) = val;
            }
        }
    }

    /*!
     * @brief Lightweight, capture-by-value view: the folded base pointer and
     * strides copied out once, so a loop body indexes local values instead of
     * reloading this->m_base/m_js/m_ks through operator(). Mirrors the amrex::Array4 captured by FIELDLOOP under USE_AMREX
     * (materialized once per tile, then indexed by value) so the same LOOP
     * macros can be hoisted out of the loop body and stay vectorizable.
     */
    struct View
    {
        T *__restrict base;        ///< already origin-folded, see cache_addressing()
        stride_t js, ks;
        inline T &operator()(int ii, int jj, int kk) const noexcept
        { return base[ii*js + jj*ks + kk]; }
    };

    struct ConstView
    {
        const T *__restrict base;
        stride_t js, ks;
        inline const T &operator()(int ii, int jj, int kk) const noexcept
        { return base[ii*js + jj*ks + kk]; }
    };

    inline View view() noexcept { return View{m_base, m_js, m_ks}; }
    inline ConstView const_view() const noexcept { return ConstView{m_base, m_js, m_ks}; }

    T *data() noexcept {return V.data();}
    const T *data() const noexcept {return V.data();}

    // whole array including ghost cells and slack; 0 for allocate=false
    std::size_t size() const noexcept {return V.size();}
    T *begin() noexcept {return V.data();}
    const T *begin() const noexcept {return V.data();}
    T *end() noexcept {return V.data()+V.size();}
    const T *end() const noexcept {return V.data()+V.size();}

protected:
    /*!
     * @brief Vertical-extent-parameterised constructor.
     *
     * The whole class is agnostic about the vertical stride — cache_addressing()
     * folds whatever @p kz is into m_ks/m_js — so a field with a different number
     * of z-planes needs nothing but a different kz. field7 passes p->kmaxF to get
     * the sigma-grid vertical-node layout (one plane more than p->kmax), which
     * reproduces the FIJK addressing in iterators3D.h exactly:
     *
     *   V[(i-imin)*jmax*kmaxF + (j-jmin)*kmaxF + k-kmin] == m_base[i*m_js + j*m_ks + k]
     *
     * @p slack is extra trailing elements, for layouts whose forward-stencil
     * macros reach past the last in-stride slot. See field7.
     * 
     * allocate=false: the field is not used in this run (V stays empty)
     */
    field_base(lexer *pp, int kz, std::size_t slack, bool allocate=true) :
        V(allocate ? static_cast<std::size_t>(pp->imax)*static_cast<std::size_t>(pp->jmax)*static_cast<std::size_t>(kz) + slack : 0, T{}),
        p(pp),
        imin(p->imin), jmin(p->jmin),
        kmin(p->kmin), kmax(kz),
        jkmax(static_cast<stride_t>(p->jmax)*kz)
    {cache_addressing();}

    /*!
     * @brief Precomputes the folded addressing used by operator() and view().
     *
     * The origin and both strides are folded into m_base so that
     *
     *   V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin] == m_base[ii*m_js + jj*m_ks + kk]
     *
     * imin/jmin/kmin are all -margin, so m_base lands inside the allocation
     * rather than before it.
     *
     * The strides are long for every payload — see stride_t.
     *
     * Must be called after anything that can reseat V's storage.
     */
    void cache_addressing() noexcept
    {
        if(V.empty())
        {
            // allocate=false: no folded base into an empty buffer (and none left
            // dangling if V was emptied after being allocated)
            m_base = nullptr;
            m_js   = 0;
            m_ks   = 0;
            return;
        }

        m_js   = jkmax;
        m_ks   = kmax;
        m_base = V.data() - static_cast<long>(imin)*jkmax
                          - static_cast<long>(jmin)*kmax - kmin;
    }

    std::vector<T> V;

    lexer *const p;

private:
    const int imin,jmin,kmin,kmax;
    const stride_t jkmax; ///< i-stride; stride_t like m_js, which it seeds

    T       *m_base = nullptr;
    stride_t m_js   = 0;
    stride_t m_ks   = 0;
};

#endif
