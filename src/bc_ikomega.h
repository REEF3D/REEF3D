/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"increment.h"
#include"roughness.h"
class fdm;
class lexer;
class field;

#ifndef BC_IKOMEGA_H_
#define BC_IKOMEGA_H_

using namespace std;

class bc_ikomega : public roughness
{
public:
	bc_ikomega(lexer*);
	virtual ~bc_ikomega();
	void bckomega_start(fdm*,lexer*,field&,field&, int);
    void bckin_matrix(fdm*,lexer*,field&,field&);
    void bcomega_matrix(fdm*,lexer*,field&,field&);
	void wall_law_kin(fdm*,lexer*,field&,field&,int,int,int,int,int,int,double);
	void wall_law_omega(fdm*,lexer*,field&,field&,int,int,int,int,int,int,double);

private:
	double uplus,ks_plus,dist,ks,ustar,u_abs,eps_star,tau;
	int ii,jj,kk;
	int count,q;
	double fac,value;
	double kappa;

};
#endif

