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

#ifndef NHFLOW_KEPSILON_FUNC_H_
#define NHFLOW_KEPSILON_FUNC_H_

#include"nhflow_rans_io.h"
#include"nhflow_kepsilon_bc.h"
#include"ghostcell.h"
#include"vrans.h"

using namespace std;

class nhflow_kepsilon_func : public nhflow_rans_io, public nhflow_kepsilon_bc
{
public:
	nhflow_kepsilon_func(lexer *, fdm_nhf*, ghostcell*);
	virtual ~nhflow_kepsilon_func();
	void isource(lexer*,fdm_nhf*) override;
	void jsource(lexer*,fdm_nhf*) override;
	void ksource(lexer*,fdm_nhf*) override;
	void kinsource(lexer*,fdm_nhf*,vrans*) override;
	void epssource(lexer*,fdm_nhf*,vrans*) override;
	void epsfsf(lexer*,fdm_nhf*,ghostcell*) override;
	void eddyvisc(lexer*,fdm_nhf*,ghostcell*,vrans*) override;
	void clearfield(lexer*,fdm_nhf*,double*) override;

	int count,q;
	double starttime;
    
private:
    double epsi;
	double dirac;
    double dxm,f;
};

#endif


