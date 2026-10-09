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


#ifndef SEDIMENT_FDM_H_
#define SEDIMENT_FDM_H_

#include"sliceint4.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"field4a.h"

using namespace std;

class sediment_mixture;

class sediment_fdm
{
public:
    sediment_fdm(lexer*);
	virtual ~sediment_fdm();
    
    slice1 P;
    slice2 Q;
    
    slice4 bedzh,bedzh0,bedch,bedsole;
    slice4 vz,dh,reduce;
    slice4 ks,ks_eff;
    slice4 ro;
    
    slice4 tau_eff,tau_crit;
    slice4 shearvel_eff,shearvel_crit;
    slice4 shields_eff, shields_crit;
    
    slice4 alpha,teta,gamma,beta,phi;
    sliceint4 active;
    
    
    sliceint4 bedk;
    slice4 slide_fh;
    
    slice4 qb,qbe,qbs;
    slice4 cbe,cb,cbn,conc;
    slice4 dryd;    // NHFLOW: suspended sediment volume per area of columns that fell dry, deposited at the next bed update
    
    slice4 waterlevel;
    slice4 guard;
    slice4 MOB,tau_i;
    
    double ws;
    double bedmax, bedmin;
    
    // grain diameter used by the bedload formulas; S20 for single-fraction runs,
    // set to d_k by sediment_mixture while the bedload of fraction k is evaluated
    double dk;
    
    // multi-fraction bed model (S51>0), nullptr otherwise
    sediment_mixture *pmix;

    int *DFBED;

};

#endif
