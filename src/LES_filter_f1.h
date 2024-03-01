/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"LES_filter.h"
#include"strain.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef LES_FILTER_F1_H_
#define LES_FILTER_F1_H_

class LES_filter_f1 : public LES_filter, public strain
{
public:
	LES_filter_f1(lexer *, fdm*);
	virtual ~LES_filter_f1();
    
	virtual void start(lexer*, fdm*, ghostcell*,field&,field&,field&, int);
    
//    int veleval(lexer*,int);


private:
    field1 ubar,etta_u1,etta_u2;
    field2 vbar,etta_v1,etta_v2;
    field3 wbar,etta_w1,etta_w2;


};

#endif


