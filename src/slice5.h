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

#ifndef SLICE5_H_
#define SLICE5_H_

#include"slice.h"
#include"increment.h"

using namespace std;

class slice5 final : public slice, increment
{
public:
	slice5(lexer* p) : slice(p) {};
	virtual ~slice5() = default;

    inline double& operator()(int ii, int jj) noexcept override final {return V[(ii-imin)*jmax + (jj-jmin)];};
    void ggcpol(lexer*) override final {};
};

#endif
