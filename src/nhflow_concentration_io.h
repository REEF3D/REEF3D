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

#ifndef NHFLOW_CONCENTRATION_IO_H_
#define NHFLOW_CONCENTRATION_IO_H_

#include"nhflow_concentration.h"
#include"increment.h"
#include"field4.h"
#include"fluid_update.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

class nhflow_concentration_io : public nhflow_concentration, increment
{
public:
    nhflow_concentration_io(lexer*,fdm_nhf*);
	virtual ~nhflow_concentration_io();

    void print_3D(lexer*, fdm_nhf*, ghostcell*, std::vector<char>&, size_t&) override final;
    void ini(lexer*, fdm_nhf*, ghostcell*, nhflow_concentration *pnhflow_concentration) override final;
    double val(int,int,int) override final;

    void name_ParaView_parallel(lexer*, ofstream&) override final;
    void name_ParaView(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView(lexer*, int*, int &) override final;
    
    double *C;

private:

	float ffn;
	double ddn;
	int n,iin;
	
	double fx(double,double,double,double,double);
	double fz(double,double,double,double,double);
};

#endif

