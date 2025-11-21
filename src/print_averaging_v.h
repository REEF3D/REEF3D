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

#ifndef PRINT_AVERAGING_V_H_
#define PRINT_AVERAGING_V_H_

#include"print_averaging.h"
#include"increment.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include<iostream>

using namespace std;

class print_averaging_v : public print_averaging, public increment
{
public:
    print_averaging_v(lexer*,fdm*,ghostcell*);
	virtual ~print_averaging_v();
    
    virtual void averaging(lexer *p, fdm *a, ghostcell *pgc, heat*);
    
    virtual void name_ParaView_parallel(lexer *p, ofstream &result);
    virtual void name_ParaView(lexer *p, std::stringstream &result, int *offset, int &n);
    virtual void offset_ParaView(lexer *p, int *offset, int &n);
    virtual void print_3D(lexer* p, fdm *a, ghostcell *pgc,  std::vector<char>&, size_t&);


private:




};

#endif
