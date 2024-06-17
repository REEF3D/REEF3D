/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef BEDLOAD_MPM_H_
#define BEDLOAD_MPM_H_

#include"bedload.h"
#include"increment.h"

class turbulence;

using namespace std;

class bedload_MPM : public bedload, public increment
{
public:

    bedload_MPM(lexer*);
	virtual ~bedload_MPM();
	virtual void start(lexer*, ghostcell*, sediment_fdm*);

private:
    double rhosed,rhowat;
    double g,d50;
    double shields,eta,visc;
    double kappa,u_plus,ks;
    double tau_eff, shearvel_eff, shields_eff;
    double tau_crit, shearvel_crit, shields_crit;

};

#endif

