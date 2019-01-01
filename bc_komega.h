/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"increment.h"
#include"roughness.h"
class fdm;
class lexer;
class field;

#ifndef BC_KOMEGA_H_
#define BC_KOMEGA_H_

using namespace std;

class bc_komega : private increment, public roughness
{
public:
	bc_komega(lexer*);
	virtual ~bc_komega();
	void bckeps_start(fdm*,lexer*,field&,field&, int);
	void wall_law_kin(fdm*,lexer*,field&,field&,int,int,int,int,int,double);
	void wall_law_omega(fdm*,lexer*,field&,field&,int,int,int,int,int,double);

private:
	double uplus,ks_plus,dist,ks,ustar,u_abs,eps_star,tau;
	int ii,jj,kk;
	double deltax,fac,value;
	const double kappa;

};
#endif
