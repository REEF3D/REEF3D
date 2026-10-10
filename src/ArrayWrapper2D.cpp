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

ArrayWrapper2D::ArrayWrapper2D(lexer *pp) : p(pp)
{
}

ArrayWrapper2D::~ArrayWrapper2D()
{
    // out-of-line: the USE_AMREX path gives this a body (deregister_imf)
}

void ArrayWrapper2D::resize(int default_value)
{
    m_data.resize(static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax), default_value);
    cache_addressing();
}

void ArrayWrapper2D::setVal(int val, bool includeGhost)
{
    if(m_data.empty())
    {
        return; // not resize()d yet: m_base is null, there is nothing to set
    }
    else if(includeGhost)
    {
        std::fill(m_data.begin(), m_data.end(), val);
    }
    else
    {
        int i,j;
        SLICELOOP4
        {
            operator()(i,j) = val;
        }
    }
}

/*!
    * @brief Precomputes the flat addressing used by operator() and operator[].
    *
    * 2D counterpart of ArrayWrapper3D::cache_addressing — see there for the
    * rationale. Here
    *
    *   m_data[(ii-imin)*jmax + (jj-jmin)]  ==  m_base[ii*m_js + jj]
    *
    * so operator() needs no lexer dereference. imin/jmin are -margin, so m_base
    * lands inside the allocation rather than before it.
*/
void ArrayWrapper2D::cache_addressing() noexcept
{
    m_js   = p->jmax;
    m_base = m_data.data() - p->imin*m_js - p->jmin;
}
