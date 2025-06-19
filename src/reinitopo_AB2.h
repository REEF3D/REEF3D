/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef REINITOPO_AB2_H_
#define REINITOPO_AB2_H_

#include"reinitopo.h"
#include"field4a.h"
#include"gradient.h"

class reinidisc;

using namespace std;

class reinitopo_AB2 : public reinitopo, gradient
{
public:
	reinitopo_AB2(lexer* p);
	virtual ~reinitopo_AB2();
	virtual void start(lexer*,fdm*,ghostcell*,field&);

	int *sizeM;



private:
	reinidisc *prdisc;
	field4a f,frk1,frk2,L,dt;
    
    void step(lexer*, fdm*);
    void time_preproc(lexer*);

	int gcval_topo,gcval_initopo,reiniter,gcval;	
};

#endif
