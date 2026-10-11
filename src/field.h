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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef FIELD_H_
#define FIELD_H_

#include "field_base.h"

#include <cassert>

class field : public field_base<double>
{
public:
    field(lexer *pp, bool allocate=true) : field_base<double>(pp,allocate) {}
    ~field() override;

    // same layout only; copies in place, so V is never reallocated and the
    // folded addressing stays valid
    void CopyFrom(const field &src)
    {
        assert(size()==src.size());
        std::ranges::copy(src, begin());
    }

protected:
    field(lexer *pp, int kz, std::size_t slack) : field_base<double>(pp, kz, slack) {}
};

#endif
