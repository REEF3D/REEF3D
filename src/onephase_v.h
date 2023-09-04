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

using namespace std;

#ifndef ONEPHASE_V_H_
#define ONEPHASE_V_H_

class onephase_v : public onephase
{
public:
    onephase_v(lexer*, fdm*, ghostcell*);
	virtual ~onephase_v();
    
	virtual void update(lexer*, fdm*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*);
    
    virtual void uvel(lexer*, fdm*, ghostcell*, field&);
    virtual void vvel(lexer*, fdm*, ghostcell*, field&);
    virtual void wvel(lexer*, fdm*, ghostcell*, field&);

};

#endif
