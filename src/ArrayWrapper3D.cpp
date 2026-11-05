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

#include "ArrayWrapper3D.h"
#include "ArrayWrapper3D_imp.h"
#include "lexer.h"
#include <algorithm>
#include <cassert>

ArrayWrapper3D::ArrayWrapper3D(lexer *pp, unsigned int loc) : p(pp), data_location(loc)
{
    assert((data_location >= 1 && data_location <= 4) || data_location == 7);
}

ArrayWrapper3D::~ArrayWrapper3D() = default;

int ArrayWrapper3D::kz() const noexcept
{
    return data_location == 7 ? p->kmaxF : p->kmax;
}

// Call after the grid metrics are final. Safe to repeat: every call refills the whole
// array (ghosts included) with default_value, so no value from an earlier size
// survives at the wrong strides.
void ArrayWrapper3D::resize(int default_value)
{
    // Single level: one flat array, sized from the same lexer metrics the IJK
    // (or, for the vertical-node layout, FIJK) macro reads.
    // grid::assign_margin sets imax/jmax/kmax/kmaxF and imin/jmin/kmin
    // together, so all are final by the time resize runs.
    //
    // The slack plane matches the imax*jmax*(kmax+2) allocation this replaces
    // in driver_makegrid_sigma.cpp: the stride is kmaxF = kmax+1, and the
    // forward-stencil macros (FIJKp3/p4) reach past the last in-stride slot in
    // the final column. See field7 for the same reasoning.
    const std::size_t plane = static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax);
    const std::size_t slack = (data_location == 7) ? plane : 0;
    m_data.assign(plane*static_cast<std::size_t>(kz()) + slack, default_value);
    cache_addressing();
}

void ArrayWrapper3D::release() noexcept
{
    std::vector<int>().swap(m_data);
    m_base = nullptr;
    m_js = 0;
    m_ks = 0;
}

/*!
 * @brief Precomputes the flat addressing used by operator() and operator[].
 *
 * See field_base::cache_addressing for the rationale (origin folding, why
 * m_base stays inside the allocation, stride width). Here
 *
 *   m_data[(i-imin)*jmax*kz + (j-jmin)*kz + (k-kmin)]  ==  m_base[i*m_js + j*m_ks + k]
 *
 * with kz = kz(): p->kmax, or p->kmaxF for the vertical-node layout (7),
 * where this is the FIJK addressing.
 *
 * Specific to this class: operator() needs no lexer dereference at all. The
 * wrapper's own p is a different pointer to the compiler than the p at the
 * call site, so nothing computed through it could be shared with neighbouring
 * accesses; the cached strides can.
 */
void ArrayWrapper3D::cache_addressing() noexcept
{
    assert(p->imin<=0 && p->jmin<=0 && p->kmin<=0); // keeps m_base inside the allocation

    m_ks   = kz();
    m_js   = static_cast<stride_t>(p->jmax) * kz();
    m_base = m_data.data() - p->imin*m_js - p->jmin*m_ks - p->kmin;
}
