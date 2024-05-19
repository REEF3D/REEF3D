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

#include"rheology.h"

using namespace std;

#ifndef RHEOLOGY_V_H_
#define RHEOLOGY_V_H_


class rheology_v : public rheology
{

public:

	rheology_v(lexer*,fdm*);
	virtual ~rheology_v();

    virtual double viscosity(lexer*,fdm*,ghostcell*);
    
    virtual void u_source(lexer*,fdm*);
    virtual void v_source(lexer*,fdm*);
    virtual void w_source(lexer*,fdm*);
    
    virtual void filltau(lexer*,fdm*,ghostcell*);

private:
    double Herschel_Bulkley(lexer*,fdm*,ghostcell*);
    
    double val;
	
};
#endif
