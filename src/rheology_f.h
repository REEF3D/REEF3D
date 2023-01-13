/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"rheology.h"
#include"strain.h"
#include"field4.h"

using namespace std;

#ifndef RHEOLOGY_F_H_
#define RHEOLOGY_F_H_

class rheology_f : public rheology, public strain
{

public:

	rheology_f(lexer*, fdm*);
	virtual ~rheology_f();

    virtual double viscosity(lexer*,fdm*,ghostcell*);
    
    virtual void u_source(lexer*,fdm*);
    virtual void v_source(lexer*,fdm*);
    virtual void w_source(lexer*,fdm*);
    
    virtual void filltau(lexer*,fdm*,ghostcell*);

private:
    field4 tau_x,tau_y,tau_z;
    
    double Herschel_Bulkley(lexer*,fdm*,ghostcell*);
    
	
	double gamma;
    double val,f,H,phival,pval;
    double tau0;
    double tau0_p,tau0_m;
    double tanphi;
    
    const double epsi;
    
    int count;
	
};
#endif
