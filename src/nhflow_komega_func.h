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

#ifndef NHFLOW_IKOMEGA_H_
#define NHFLOW_IKOMEGA_H_

#include"nhflow_rans_io.h"
#include"nhflow_komega_bc.h"
#include"ghostcell.h"
#include"vrans.h"

using namespace std;

class nhflow_komega_func : public nhflow_rans_io, public nhflow_komega_bc
{
public:
	nhflow_komega_func(lexer *, fdm_nhf*, ghostcell*);
	virtual ~nhflow_komega_func();
	virtual void isource(lexer*,fdm_nhf*);
	virtual void jsource(lexer*,fdm_nhf*);
	virtual void ksource(lexer*,fdm_nhf*);
	virtual void kinsource(lexer*,fdm_nhf*,vrans*);
	virtual void epssource(lexer*,fdm_nhf*,vrans*);
	virtual void epsfsf(lexer*,fdm_nhf*,ghostcell*);
	virtual void eddyvisc(lexer*,fdm_nhf*,ghostcell*,vrans*);
	virtual void clearfield(lexer*,fdm_nhf*,double*);

	int count,q;
	double starttime;
    
private:
    double epsi;
	double dirac;
};

#endif


