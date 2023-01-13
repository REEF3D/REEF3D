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

#include"turbulence.h"
#include"increment.h"

using namespace std;

#ifndef KEPSILON_VOID_H_
#define KEPSILON_VOID_H_

class kepsilon_void : public turbulence, public increment
{

public:
	kepsilon_void(lexer *,fdm*,ghostcell*);
	virtual ~kepsilon_void();

	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*);
	virtual void ktimesave(lexer*, fdm*, ghostcell*);
	virtual void etimesave(lexer*, fdm*, ghostcell*);

	virtual void isource(lexer*, fdm*);
	virtual void jsource(lexer*, fdm*);
	virtual void ksource(lexer*, fdm*);

	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void ini(lexer*, fdm*, ghostcell*);
    virtual double kinval(int,int,int);
    virtual double epsval(int,int,int);
	virtual double ccipol_kinval(lexer*,ghostcell*,double,double,double);
	virtual double ccipol_epsval(lexer*,ghostcell*,double,double,double);
    virtual double ccipol_a_kinval(lexer*,ghostcell*,double,double,double);
	virtual double ccipol_a_epsval(lexer*,ghostcell*,double,double,double);
    virtual void kinget(int,int,int,double);
    virtual void epsget(int,int,int,double);
	virtual void gcupdate(lexer*, fdm*, ghostcell*);
	
    virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
};

#endif
