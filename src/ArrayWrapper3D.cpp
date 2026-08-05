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

ArrayWrapper3D::ArrayWrapper3D(lexer *pp, unsigned int _data_location) : p(pp), data_location(_data_location)
{
}

ArrayWrapper3D::~ArrayWrapper3D()
{
    // out-of-line: the USE_AMREX path gives this a body (deregister_imf)
}

void ArrayWrapper3D::resize(int default_value)
{
    // Single level: one flat array, sized from the same lexer metrics the IJK
    // macro reads. grid::assign_margin sets imax/jmax/kmax and imin/jmin/kmin
    // together, so both are final by the time resize runs.
    m_data.resize(static_cast<std::size_t>(p->imax)*p->jmax*p->kmax, default_value);
    cache_addressing();
}

void ArrayWrapper3D::setVal(int val, bool includeGhost)
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
        // 1/2/3 included: same layout as 4 without AMReX (see the header)
        int i,j,k;
        LOOP
        {
            operator()(i,j,k) = val;
        }
    }
}

/*!
    * @brief Precomputes the flat addressing used by operator() and operator[].
    *
    * The origin and both strides are folded into m_base, so that
    *
    *   m_data[(i-imin)*jmax*kz + (j-jmin)*kz + (k-kmin)]  ==  m_base[i*m_js + j*m_ks + k]
    *
    * and operator() needs no lexer dereference at all. That is what lets the
    * index be shared with neighbouring accesses instead of being rebuilt from
    * p on every one — the wrapper's own p is a different pointer to the
    * compiler than the p at the call site, so nothing computed through it can
    * be reused.
    *
    * m_base is not a pointer formed before the array: imin/jmin/kmin are all
    * -margin, so it lands at data.data() + margin*(m_js + m_ks + 1).
    *
    * The strides are long deliberately. An int store through the reference
    * returned by operator() may alias an int member under TBAA, which forces a
    * reload of the strides on every iteration of a writing loop; as longs they
    * stay in registers and the address strength-reduces to an induction
    * variable.
*/
void ArrayWrapper3D::cache_addressing() noexcept
{
    m_ks   = p->kmax;
    m_js   = static_cast<long>(p->jmax) * p->kmax;
    m_base = m_data.data() - p->imin*m_js - p->jmin*m_ks - p->kmin;
}
