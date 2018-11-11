/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"rans_io.h"
#include"bc_ikomega.h"
#include"ghostcell.h"

class multiphase;
class vrans;

using namespace std;

#ifndef IKOMEGA_H_
#define IKOMEGA_H_

class ikomega : public rans_io, public bc_ikomega
{
public:
	ikomega(lexer *, fdm*, ghostcell*, multiphase*);
	virtual ~ikomega();
	virtual void isource(lexer*,fdm*);
	virtual void jsource(lexer*,fdm*);
	virtual void ksource(lexer*,fdm*);
	virtual void kinsource(lexer*,fdm*);
	virtual void epssource(lexer*,fdm*);
	virtual void epsfsf(lexer*,fdm*,ghostcell*);
	virtual void eddyvisc(lexer*,fdm*,ghostcell*);
	virtual void clearfield(lexer*,fdm*,field&);

	int count,q;
	double starttime;
	
	multiphase *pmp;	
    vrans *pvrans;
};

#endif


