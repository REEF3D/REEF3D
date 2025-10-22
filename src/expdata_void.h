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

#ifndef EXPDATA_VOID_H_
#define EXPDATA_VOID_H_

#include"expdata.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

class expdata_void : public expdata
{
public:
	expdata_void(lexer*, fdm*, ghostcell*);
	virtual ~expdata_void();
	virtual void start(lexer*, fdm*, ghostcell*);
	
	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_ParaView(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_ParaView(lexer*, int*, int &);
};

#endif
