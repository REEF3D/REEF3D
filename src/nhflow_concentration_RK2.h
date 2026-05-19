/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef NHFLOW_CONCENTRATION_RK2_H_
#define NHFLOW_CONCENTRATION_RK2_H_

#include"nhflow_concentration_io.h"

using namespace std;

class nhflow_concentration_RK2 final : public nhflow_concentration_io
{
public:
    nhflow_concentration_RK2(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_concentration_RK2();

	void start(lexer*, fdm_nhf*, convection*, diffusion*, turbulence*, solver*, ghostcell*, ioflow*) override final;
	void ctimesave(lexer*, fdm_nhf*) override final;


private:
    void clearrhs(lexer*,fdm_nhf*,ghostcell*);
    
	int gcval_concentration;
	double starttime, endtime;
};

#endif
