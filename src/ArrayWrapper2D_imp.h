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

#ifndef ARRAYWRAPPER2D_IMP_H_
#define ARRAYWRAPPER2D_IMP_H_

// Inline hot-path methods of ArrayWrapper2D. Included at the bottom of
// lexer.h, after the lexer class is complete (same pattern as
// ArrayWrapper3D_imp.h). The bodies below do not touch p (the addressing is
// cached), so they only need the declaration, not the complete lexer.

#include "ArrayWrapper2D.h"
#include "lexer.h"

// Folded addressing, see field_base::cache_addressing. These touch no lexer
// member. Equivalent to data[(ii-imin)*jmax + (jj-jmin)].
inline int &ArrayWrapper2D::operator()(int ii, int jj) noexcept
{
    return m_base[ii*m_js + jj];
}

inline const int &ArrayWrapper2D::operator()(int ii, int jj) const noexcept
{
    return m_base[ii*m_js + jj];
}

// See ArrayWrapper3D::operator[].
inline int &ArrayWrapper2D::operator[](std::size_t index) noexcept
{
    return m_data[index];
}

inline const int &ArrayWrapper2D::operator[](std::size_t index) const noexcept
{
    return m_data[index];
}

#endif
