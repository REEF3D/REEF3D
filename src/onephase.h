/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

class lexer;
class fdm;
class ghostcell;
class ioflow;
class field;

using namespace std;

#ifndef ONEPHASE_H_
#define ONEPHASE_H_

class onephase
{
public:
	virtual void update(lexer*, fdm*, ghostcell*, ioflow*)=0;
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*)=0;
    virtual void uvel(lexer*, fdm*, ghostcell*, field&)=0;
    virtual void vvel(lexer*, fdm*, ghostcell*, field&)=0;
    virtual void wvel(lexer*, fdm*, ghostcell*, field&)=0;

};

#endif
