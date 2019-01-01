/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"strain.h"
#include"rans_io.h"
#include"bc_komega.h"

class ghostcell;
class vrans;

using namespace std;

#ifndef KOMEGA_H_
#define KOMEGA_H_

class komega : public rans_io, public bc_komega
{
public:
	komega(lexer *, fdm*, ghostcell *pgc);
	virtual ~komega();
	virtual void isource(lexer*,fdm*);
	virtual void jsource(lexer*,fdm*);
	virtual void ksource(lexer*,fdm*);
	virtual void kinsource(lexer*,fdm*,ghostcell*,field&,field&);
	virtual void epssource(lexer*,fdm*,ghostcell*,field&,field&);
	virtual void epsfsf(lexer*,fdm*,ghostcell*,field&,field&);
	virtual void eddyvisc(lexer*,fdm*,field&,field&, ghostcell*);
	virtual void clearfield(lexer*,fdm*,field&);
    virtual void rhs(lexer*,fdm*,ghostcell*);
	
	double starttime;
    
    vrans *pvrans;
};

#endif

