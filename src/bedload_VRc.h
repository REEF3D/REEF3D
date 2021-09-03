/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"bedload.h"
#include"bedshear.h"
#include"bedload_noneq.h"

class turbulence;

using namespace std;

#ifndef BEDLOAD_VRC_H_
#define BEDLOAD_VRC_H_

class bedload_VRc : public bedload, public bedshear, public bedload_noneq
{
public:

    bedload_VRc(lexer*,turbulence*);
    virtual ~bedload_VRc();

	virtual void start(lexer*, fdm*, ghostcell*, sediment_fdm*);

private:
    const double epsi;
    double rhosed,rhowat,Rstar,Ds;
    double g,d50;
    double visc;
    double kappa,u_plus,ks;
    double tau_eff, shearvel_eff, shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
};

#endif


