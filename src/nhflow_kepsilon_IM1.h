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

#ifndef NHFLOW_KEPSILON_IM1_H_
#define NHFLOW_KEPSILON_IM1_H_

#include"nhflow_kepsilon_func.h"

using namespace std;

class nhflow_kepsilon_IM1 : public nhflow_kepsilon_func
{
public:
	nhflow_kepsilon_IM1(lexer *, fdm_nhf*, ghostcell*);
	virtual ~nhflow_kepsilon_IM1();
	void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, vrans*) override;
	void ktimesave(lexer*, fdm_nhf*, ghostcell*) override;
	void etimesave(lexer*, fdm_nhf*, ghostcell*) override;
	void timesource(lexer*,fdm_nhf*,double*);
    void kinupdate(lexer*, fdm_nhf*, ghostcell*);
	void clearrhs(lexer*,fdm_nhf*);

	double  *KN,*EN;


private:
    int gcval_kin, gcval_eps;
    int count,q;
    double aii;
};

#endif

