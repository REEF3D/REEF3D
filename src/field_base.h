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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef FIELD_BASE_H_
#define FIELD_BASE_H_

#include "lexer.h"

template<typename T>
class field_base
{
public:
    field_base(lexer* p) : imin(p->imin), imax(p->imax), jmin(p->jmin), jkmax(p->jmax*p->kmax), kmin(p->kmin), kmax(p->kmax)
    {
        V = new T[imax*jkmax] {};
    }
    virtual ~field_base()
    {
        delete [] V;
        V = nullptr;
    }

    field_base(const field_base&) = delete;
    field_base& operator=(const field_base&) = delete;
    field_base(field_base&&) = delete;
    field_base& operator=(field_base&&) = delete;

    inline T& operator()(int ii, int jj, int kk) noexcept {return V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin];};
    inline T& operator[](int n) noexcept {return V[n];};

	T *V;

private:
    const int imin,imax,jkmax,jmin,kmin,kmax;
};

#endif
