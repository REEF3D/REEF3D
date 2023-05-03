/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"bedshear_reduction.h"
#include"bedslope.h"

class lexer;
class ghostcell;
class sediment_fdm;

using namespace std;

#ifndef REDUCTION_DEYANA_H_
#define REDUCTION_DEYANA_H_

class reduction_deyana :  public bedshear_reduction, public bedslope
{
public:
    reduction_deyana(lexer*);
    virtual ~reduction_deyana();

	virtual void start(lexer*,ghostcell*,sediment_fdm*);

private:

    double u_abs,u_plus,dist;
    double uvel, vvel;
    double tetaval, alphaval, phival, tanphi, beta;
    double eta;
};

#endif


