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
Author: Ahmet Soydan
--------------------------------------------------------------------*/


#ifndef DENSITY_CAV_H_
#define DENSITY_CAV_H_

#include"density.h"
#include"increment.h"

class fdm;
class lexer;

// face density of the liquid-vapour mixture of the cavitation model (W 60), from the vapour
// fraction of the two cells; the liquid density W 1 without an active model
class density_cav final : public density, virtual public increment
{
public:
    density_cav(lexer*);
    virtual ~density_cav();

    double roface(lexer*,fdm*,int,int,int) override final;
};

#endif
