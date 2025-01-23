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

#include"nhflow_les_io.h"

using namespace std;

#ifndef NHFLOW_LES_SMAGORINSKY_H_
#define NHFLOW_LES_SMAGORINSKY_H_

class nhflow_LES_Smagorinsky : public nhflow_les_io
{
public:
	nhflow_LES_Smagorinsky(lexer *, fdm_nhf*, ghostcell*);
	virtual ~nhflow_LES_Smagorinsky();
	virtual void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, vrans*);
	virtual void ktimesave(lexer*, fdm_nhf*, ghostcell*);
	virtual void etimesave(lexer*, fdm_nhf*, ghostcell*);
	void timesource(lexer*,fdm_nhf*,double*);
    void kinupdate(lexer*, fdm_nhf*, ghostcell*);
	void clearrhs(lexer*,fdm_nhf*);

private:
    int gcval_kin, gcval_eps;
    int count,q;
    double c_sgs;
    double aii;
};

#endif

