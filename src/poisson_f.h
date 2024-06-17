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

#ifndef POISSON_F_H_
#define POISSON_F_H_

#include"poisson.h"
#include"increment.h"

class heat;
class concentration;
class density;

using namespace std;

class poisson_f : public poisson, public increment
{

public:

	poisson_f (lexer*, heat*&, concentration*&);
	virtual ~poisson_f();

	virtual void start(lexer *,fdm*,field&);

private:

	int count,n,q;
    
    density *pd;
};


#endif



