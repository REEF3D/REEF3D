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

#ifndef KEPSILON_VOID_H_
#define KEPSILON_VOID_H_

#include"turbulence.h"
#include"increment.h"

using namespace std;

class kepsilon_void : public turbulence, public increment
{

public:
	kepsilon_void(lexer *,fdm*,ghostcell*);
	virtual ~kepsilon_void();

	void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*) override final;
	void ktimesave(lexer*, fdm*, ghostcell*) override final;
	void etimesave(lexer*, fdm*, ghostcell*) override final;

	void isource(lexer*, fdm*) override final;
	void jsource(lexer*, fdm*) override final;
	void ksource(lexer*, fdm*) override final;

	void print_3D(lexer*, fdm*, ghostcell*, std::vector<char>&, size_t&) override final;
    void ini(lexer*, fdm*, ghostcell*) override final;
    double kinval(int,int,int) override final;
    double epsval(int,int,int) override final;
	double ccipol_kinval(lexer*,ghostcell*,double,double,double) override final;
	double ccipol_epsval(lexer*,ghostcell*,double,double,double) override final;
    double ccipol_a_kinval(lexer*,ghostcell*,double,double,double) override final;
	double ccipol_a_epsval(lexer*,ghostcell*,double,double,double) override final;
    void kinget(int,int,int,double) override final;
    void epsget(int,int,int,double) override final;
	void gcupdate(lexer*, fdm*, ghostcell*) override final;
	
    void name_ParaView_parallel(lexer*, ofstream&) override final;
    void name_ParaView(lexer*, std::stringstream&, int*, int &) override final;
    void offset_ParaView(lexer*, int*, int &) override final;
};

#endif
