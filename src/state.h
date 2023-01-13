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

#include"increment.h"


class lexer;
class fdm;
class ghostcell;
class turbulence;
class sediment;

using namespace std;

#ifndef STATE_H_
#define STATE_H_

class state : public increment
{

public:
	state(lexer*,fdm*,ghostcell*);
	virtual ~state();
	virtual void write(lexer*,fdm*,ghostcell*,turbulence*,sediment*);
    virtual void read(lexer*,fdm*,ghostcell*,turbulence*,sediment*);
	
private:
    virtual void filename(lexer*,fdm*,ghostcell*,int);

    char name[200];
    float ffn;
	int iin;
	double ddn;
	int printcount;
    
    
};

#endif



