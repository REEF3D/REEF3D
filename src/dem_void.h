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

#ifndef DEM_VOID_H_
#define DEM_VOID_H_

#include"dem.h"

class dem_void final : public dem
{
public:
    dem_void() = default;
    virtual ~dem_void() = default;

    void start_cfd(lexer*, fdm*, ghostcell*) override final {}
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*) override final {}
    void forcing_cfd(lexer*, fdm*, ghostcell*, int, double, field&, field&, field&, bool) override final {}
    void forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, int, double, double*, double*, double*, slice&, bool, bool) override final {}
};

#endif
