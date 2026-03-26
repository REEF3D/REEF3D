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

#ifndef MULTIPHASE_V_H_
#define MULTIPHASE_V_H_

class fdm;
class lexer;
class convection;
class solver;
class ghostcell;
class ioflow;
class reini;
class field;

#include"multiphase.h"
#include<fstream>

using namespace std;

class multiphase_v final : public multiphase
{
public:
	multiphase_v();
	virtual ~multiphase_v();
	void start(lexer*,fdm*,ghostcell*,convection*,solver*,ioflow*,reini*,particle_corr*) override final;
	void ini(lexer*,fdm*,ghostcell*,ioflow*,convection*,solver*) override final;
	void update(lexer*,fdm*,ghostcell*) override final;
	
	void print_3D(lexer*, fdm*, ghostcell*, std::vector<char>&, size_t&) override final;
	void print_file(lexer*, fdm*, ghostcell*) override final;
    double ls1val(int,int,int) override final;
    double ls2val(int,int,int) override final;
	double ccipol_ls1val(lexer*,ghostcell*,double,double,double) override final;
	double ccipol_ls2val(lexer*,ghostcell*,double,double,double) override final;
    void ls1get(int,int,int,double) override final;
    void ls2get(int,int,int,double) override final;

    void name_ParaView_parallel(lexer*, ofstream&) override final;
    void name_ParaView(lexer*, std::stringstream&, int*, int &) override final;
    void offset_ParaView(lexer*, int*, int &) override final;
};

#endif
