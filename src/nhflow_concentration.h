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

#ifndef NHFLOW_CONCENTRATION_H_
#define NHFLOW_CONCENTRATION_H_

class fdm_nhf;
class lexer;
class convection;
class diffusion;
class solver;
class ghostcell;
class ioflow;
class turbulence;
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

using namespace std;

class nhflow_concentration
{
public:

	virtual void start(lexer*, fdm_nhf*, convection*, diffusion*, turbulence*, solver*, ghostcell*, ioflow*)=0;
	virtual void ini(lexer*, fdm_nhf*, ghostcell*, nhflow_concentration*)=0;
	virtual void ctimesave(lexer*, fdm_nhf*)=0;

	virtual void print_3D(lexer*, fdm_nhf*, ghostcell*, std::vector<char>&, size_t&)=0;
	virtual double val(int,int,int)=0;

    virtual void name_ParaView_parallel(lexer*, ofstream&)=0;
    virtual void name_ParaView(lexer*, ostream&, int*, int &)=0;
    virtual void offset_ParaView(lexer*, int*, int &)=0;
};

#endif
