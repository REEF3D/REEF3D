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

#include"vorticity.h"
#include<fstream>

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef VORTICITY_VOID_H_
#define VORTICITY_VOID_H_

class vorticity_void : public vorticity
{
public:
    vorticity_void(lexer*,fdm*);
	virtual ~vorticity_void();

    virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);

    virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);

};

#endif




