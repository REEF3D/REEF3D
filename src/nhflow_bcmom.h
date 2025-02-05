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

#include"surftens.h"
#include"roughness.h"
class lexer;
class fdm_nhf;
class ghostcell;
class field;
class turbulence;

#ifndef BCMOM_NHFLOW_H_
#define BCMOM_NHFLOW_H_

using namespace std;

class nhflow_bcmom : public roughness
{
public:
	nhflow_bcmom(lexer*);
	virtual ~nhflow_bcmom();
	virtual void nhflow_bcmom_start(fdm*,lexer*,ghostcell*,turbulence*,field&, int);
	void roughness_u(lexer*, fdm_nhf*, double*, double*, slice&);
	void roughness_v(lexer*, fdm_nhf*, double*, double*, slice&);
	void roughness_w(lexer*, fdm_nhf*, double*, double*, slice&);

private:
	const double kappa;
	double uplus,ks_plus,dist,ks,ustar;
	int ii,jj,kk;
	double value;
	int gcval_phi, bckin, wallfunc_type;
};
#endif
