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

#include "ArrayWrapper2D.h"
#include "ArrayWrapper2D_imp.h"
#include "lexer.h"
#include <algorithm>
#include <cassert>

ArrayWrapper2D::ArrayWrapper2D(lexer *pp) : p(pp)
{
}

ArrayWrapper2D::~ArrayWrapper2D() = default;

// Call after the grid metrics are final. Safe to repeat: every call refills the
// whole array (ghosts included) with default_value, so no value from an earlier
// size survives at the wrong strides.
void ArrayWrapper2D::resize(int default_value)
{
    m_data.assign(static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax), default_value);
    cache_addressing();
}

void ArrayWrapper2D::release() noexcept
{
    std::vector<int>().swap(m_data);
    m_base = nullptr;
    m_js = 0;
}

/*!
 * @brief Precomputes the flat addressing used by operator() and operator[].
 *
 * See field_base::cache_addressing for the rationale (origin folding, why
 * m_base stays inside the allocation, stride width). Here
 *
 *   m_data[(ii-imin)*jmax + (jj-jmin)]  ==  m_base[ii*m_js + jj]
 */
void ArrayWrapper2D::cache_addressing() noexcept
{
    assert(p->imin<=0 && p->jmin<=0); // keeps m_base inside the allocation

    m_js   = p->jmax;
    m_base = m_data.data() - p->imin*m_js - p->jmin;
}
