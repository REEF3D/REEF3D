/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"onephase.h"
#include"increment.h"

using namespace std;

#ifndef ONEPHASE_F_H_
#define ONEPHASE_F_H_

class onephase_f : public onephase, public increment
{
public:
    onephase_f(lexer*, fdm*, ghostcell*);
	virtual ~onephase_f();
    
	virtual void update(lexer*, fdm*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*);
    
private: 
    void fsf_update(lexer*,fdm*,ghostcell*);

};

#endif
