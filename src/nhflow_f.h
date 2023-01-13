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

#include"nhflow.h"
#include"increment.h"

using namespace std;

#ifndef NHFLOW_F_H_
#define NHFLOW_F_H_

class nhflow_f : public nhflow, public increment
{
public:    
    nhflow_f(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_f();

    virtual void ini(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    
    virtual void kinematic_fsf(lexer*, fdm_nhf*, double*, double*, double*, slice&, slice&, double);

private:
    int q,margin;
        

};

#endif
