/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"reduction.h"
#include"norm_vec.h"

class lexer;
class fdm;
class ghostcell;
class ddweno_f_nug;

using namespace std;

#ifndef BEDSLOPE_H_
#define BEDSLOPE_H_

class bedslope :  public norm_vec
{
public:
    bedslope(lexer*);
    virtual ~bedslope();

	virtual void slope(lexer*, fdm*,ghostcell*,double&,double&,double&,double&);

private:

    double u_abs,u_plus,dist;
    double uvel, vvel;
    double midphi, delta,beta;
    double alpha0, teta0;
    
    ddweno_f_nug *pdx;
};

#endif


