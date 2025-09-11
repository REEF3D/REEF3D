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

#ifndef IKEPSILON_H_
#define IKEPSILON_H_

#include"rans_io.h"
#include"kepsilon_bc.h"
#include"vrans.h"

using namespace std;

class kepsilon_func : public rans_io, public kepsilon_bc
{
public:
	kepsilon_func(lexer*,fdm*,ghostcell*);
	virtual ~kepsilon_func();
	void isource(lexer*,fdm*) override;
	void jsource(lexer*,fdm*) override;
	void ksource(lexer*,fdm*) override;
	void kinsource(lexer*,fdm*,vrans*);
	void epssource(lexer*,fdm*,vrans*);
	void epsfsf(lexer*,fdm*,ghostcell*);
	void eddyvisc(fdm*,lexer*,ghostcell*,vrans*);
	void clearfield(lexer*,fdm*,field&);

	int count,q;
	double starttime;
};

#endif



