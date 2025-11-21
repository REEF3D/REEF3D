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

#ifndef BCMOM_SFLOW_H_
#define BCMOM_SFLOW_H_

#include"surftens.h"
#include"roughness.h"
class lexer;
class fdm2D;
class ghostcell;
class field;
class turbulence;

using namespace std;

class sflow_bcmom : public roughness
{
public:
	sflow_bcmom(lexer*);
	virtual ~sflow_bcmom();
	virtual void sflow_bcmom_start(fdm*,lexer*,ghostcell*,turbulence*,field&, int);
	void roughness_u(lexer*, fdm2D*, slice&, slice&, slice&);
	void roughness_v(lexer*, fdm2D*, slice&, slice&, slice&);
	void roughness_w(lexer*, fdm2D*, slice&, slice&, slice&);

private:
	const double kappa;
	double uplus,ks_plus,dist,ks,ustar;
	int ii,jj,kk;
	double value;
	int gcval_phi, bckin, wallfunc_type;
};
#endif
