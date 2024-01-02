/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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


#include"sliceint4.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"field4a.h"

using namespace std;

#ifndef SEDIMENT_FDM_H_
#define SEDIMENT_FDM_H_

class sediment_fdm
{
public:
    sediment_fdm(lexer*);
	virtual ~sediment_fdm();
    
    slice1 P;
    slice2 Q;
    
    slice4 bedzh,bedzh0,bedch,bedsole;
    slice4 vz,dh,reduce;
    slice4 ks;
    
    slice4 tau_eff,tau_crit;
    slice4 shearvel_eff,shearvel_crit;
    slice4 shields_eff, shields_crit;
    
    slice4 alpha,teta,gamma,beta,phi;
    slice4 active;
    
    
    sliceint4 bedk;
    slice4 slideflag;
    
    slice4 qb,qbe;
    slice4 cbe,cb,cbn,conc;
    
    slice4 waterlevel;
    slice4 guard;
    
    double ws;

};

#endif