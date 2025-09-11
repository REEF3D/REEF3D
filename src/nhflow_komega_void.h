/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef NHFLOW_KOMEGA_VOID_H_
#define NHFLOW_KOMEGA_VOID_H_

#include"nhflow_turbulence.h"
#include"increment.h"

using namespace std;

class nhflow_komega_func_void : public nhflow_turbulence, public increment
{

public:
	nhflow_komega_func_void(lexer *,fdm_nhf*,ghostcell*);
	virtual ~nhflow_komega_func_void();

	void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, vrans*) override;
	void ktimesave(lexer*, fdm_nhf*, ghostcell*) override;
	void etimesave(lexer*, fdm_nhf*, ghostcell*) override;

	void isource(lexer*, fdm_nhf*) override;
	void jsource(lexer*, fdm_nhf*) override;
	void ksource(lexer*, fdm_nhf*) override;

    void print_2D(lexer*, fdm_nhf*, ghostcell*,ofstream&,int) override {};
	void print_3D(lexer*, fdm_nhf*, ghostcell*,ofstream&) override;
    void ini(lexer*, fdm_nhf*, ghostcell*) override;
    double kinval(int,int,int) override;
    double epsval(int,int,int) override;
	double ccipol_kinval(lexer*,ghostcell*,double,double,double) override;
	double ccipol_epsval(lexer*,ghostcell*,double,double,double) override;
    double ccipol_a_kinval(lexer*,ghostcell*,double,double,double) override;
	double ccipol_a_epsval(lexer*,ghostcell*,double,double,double) override;
    void kinget(int,int,int,double) override;
    void epsget(int,int,int,double) override;
	void gcupdate(lexer*, fdm_nhf*, ghostcell*) override;
	
    void name_pvtu(lexer*, fdm_nhf*, ghostcell*,ofstream&) override;
    void name_vtu(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtu(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &) override;
    
    void name_pvtp(lexer*, fdm_nhf*, ghostcell*,ofstream&) override {};
    void name_vtp(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtp(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &) override {};
};

#endif
