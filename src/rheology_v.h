/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef RHEOLOGY_V_H_
#define RHEOLOGY_V_H_

#include"rheology.h"

class rheology_v : public rheology
{

public:

<<<<<<< HEAD
    rheology_v();
    virtual ~rheology_v();
=======
	rheology_v(lexer*,fdm*);
	virtual ~rheology_v();
>>>>>>> parent of 516fad2a7 (Replaced \t with 4 spaces)

    double viscosity(lexer*,fdm*,ghostcell*) override;
    
    void u_source(lexer*,fdm*) override;
    void v_source(lexer*,fdm*) override;
    void w_source(lexer*,fdm*) override;
    
<<<<<<< HEAD
    void filltau(lexer*,fdm*,ghostcell*) override;
    
=======
    virtual void filltau(lexer*,fdm*,ghostcell*);

private:
    double Herschel_Bulkley(lexer*,fdm*,ghostcell*);
    
    double val;
	
>>>>>>> parent of 516fad2a7 (Replaced \t with 4 spaces)
};
#endif
