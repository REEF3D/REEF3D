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

#ifndef ARRAYWRAPPER3D_IMP_H_
#define ARRAYWRAPPER3D_IMP_H_

#include "ArrayWrapper3D.h"

inline int &ArrayWrapper3D::operator()(int i, int j, int k) noexcept
{
    // Folded addressing, see field_base::cache_addressing. Touches no lexer
    // member. Equivalent to data[IJK].
    return m_base[i*m_js + j*m_ks + k];
}

inline const int &ArrayWrapper3D::operator()(int i, int j, int k) const noexcept
{
    // Folded addressing, see field_base::cache_addressing. Touches no lexer
    // member. Equivalent to data[IJK].
    return m_base[i*m_js + j*m_ks + k];
}

// ArrayWrapper3D::operator[] still exists to provide the same interface as legacy code expects.
// Once the legacy code is update to use operator() instead, this operator[] can be removed.
inline int &ArrayWrapper3D::operator[](std::size_t index) noexcept
{
    return m_data[index];
}

inline const int &ArrayWrapper3D::operator[](std::size_t index) const noexcept
{
    return m_data[index];
}

#endif
