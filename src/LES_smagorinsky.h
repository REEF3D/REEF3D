/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"LES.h"
#include"field4.h"

class LES_filter;

using namespace std;

#ifndef LES_SMAGORINSKY_H_
#define LES_SMAGORINSKY_H_

class LES_smagorinsky : public LES
{
public:
	LES_smagorinsky(lexer *, fdm*);
	virtual ~LES_smagorinsky();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*);
	virtual void ktimesave(lexer*, fdm*, ghostcell*);
	virtual void etimesave(lexer*, fdm*, ghostcell*);

private:
	int gcval_sgs;
	double c_sgs;
	int gcval_u1, gcval_v1, gcval_w1;
	int gcval_u2, gcval_v2, gcval_w2;
    
    LES_filter *pfilter;

};

#endif


