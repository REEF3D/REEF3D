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

#include"norm_vec.h"
#include"slice4.h"

class lexer;
class fdm;
class ghostcell;
class sediment_fdm;
class turbulence;
class sliceint;

using namespace std;

#ifndef BEDSHEAR_H_
#define BEDSHEAR_H_

class bedshear :  public norm_vec
{
public:
    bedshear(lexer*,turbulence*);
    virtual ~bedshear();

	virtual void taubed(lexer*, fdm*,ghostcell*,sediment_fdm*);
	virtual void taucritbed(lexer*, fdm*,ghostcell*,sediment_fdm*);
    
    virtual void taubed(lexer*, fdm*, ghostcell*, double&);
    virtual void taucritbed(lexer*, fdm*, ghostcell*, double&);
    
    virtual void taubed(lexer*, fdm2D*,ghostcell*,sediment_fdm*);
    virtual void taucritbed(lexer*, fdm2D*,ghostcell*,sediment_fdm*);

	const double ks,kappa;
    
    slice4 taueff_loc,taucrit_loc;

private:
    turbulence *pturb;
    double tau,tauc;
    double u_abs,u_plus,dist;
    double nx,ny,nz,norm;
    double xip,yip,zip;
    double uvel,vvel,wvel;
	double kinval;
    double scale;
};

#endif
