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

#ifndef SLICE2_H_
#define SLICE2_H_

#include"slice.h"
#include"increment.h"

class slice2 final : public slice, increment
{
public:
    slice2(lexer*);
    virtual ~slice2();

    inline double& operator()(int, int) override final;
    void ggcpol(lexer*) override final;

private:
    void fieldgcalloc(lexer*);

    int iter;
    int gcfeldsize;

    int di,dj;

    lexer *pp;

    double ***gcfeld;
};

#endif
