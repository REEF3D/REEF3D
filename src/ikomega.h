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

#include"rans_io.h"
#include"bc_ikomega.h"
#include"ghostcell.h"
#include"vrans.h"

using namespace std;

#ifndef IKOMEGA_H_
#define IKOMEGA_H_

class ikomega : public rans_io, public bc_ikomega
{
public:
	ikomega(lexer *, fdm*, ghostcell*);
	virtual ~ikomega();
	virtual void isource(lexer*,fdm*);
	virtual void jsource(lexer*,fdm*);
	virtual void ksource(lexer*,fdm*);
	virtual void kinsource(lexer*,fdm*,vrans*);
	virtual void epssource(lexer*,fdm*,vrans*,field&);
	virtual void epsfsf(lexer*,fdm*,ghostcell*);
	virtual void eddyvisc(lexer*,fdm*,ghostcell*,vrans*);
	virtual void clearfield(lexer*,fdm*,field&);

	int count,q;
	double starttime;
    
private:
    double epsi;
	double dirac;
};

#endif


