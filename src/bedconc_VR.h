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

#ifndef BEDCONC_VR_H_
#define BEDCONC_VR_H_

#include"bedconc.h"
#include"increment.h"

using namespace std;

class bedconc_VR : public bedconc, public increment
{
public:
	bedconc_VR(lexer*);
	virtual ~bedconc_VR();
	void start(lexer*,ghostcell*,sediment_fdm*);

private:
	int ii,jj,kk;
	
	double d50,ks,shields,kappa;
	double Rstar, g, visc;
	double rhosed,rhowat;
	double Ti,Ds;
    double adist,zdist,deltab;

};
#endif



