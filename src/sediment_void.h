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

#include"sediment.h"

using namespace std;

#ifndef SEDIMENT_VOID_H_
#define SEDIMENT_VOID_H_

class sediment_void : public sediment
{
public:
    sediment_void();
	virtual ~sediment_void();

	virtual void start(lexer*, fdm*, convection*, ghostcell*, ioflow*, topo*, reinitopo*, suspended*, bedload*);
	virtual void update(lexer*,fdm*,ghostcell*,ioflow*);
    virtual void relax(lexer*,fdm*,ghostcell*);
    virtual void ini(lexer*,fdm*,ghostcell*);
	virtual double bedshear_point(lexer*,fdm*,ghostcell*);
    
	virtual void print_3D_bedshear(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedshear(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu_bedshear(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedshear(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_3D_parameters(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameters(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu_parameters(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameters(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
};

#endif
