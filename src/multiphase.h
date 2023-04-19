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

class fdm;
class lexer;
class convection;
class solver;
class ghostcell;
class ioflow;
class reini;
class printer;
class field;
class particle_corr;

#include<fstream>

using namespace std;

#ifndef MULTIPHASE_H_
#define MULTIPHASE_H_

class multiphase
{
public:

	virtual void start(lexer*,fdm*,ghostcell*,convection*,solver*,ioflow*,reini*,particle_corr*,printer*)=0;
	virtual void ini(lexer*,fdm*,ghostcell*,ioflow*,printer*,convection*,solver*)=0;
	virtual void update(lexer*,fdm*,ghostcell*)=0;
	
	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&)=0;
	virtual void print_file(lexer*, fdm*, ghostcell*)=0;
	virtual void nodefill(lexer*,fdm*,ghostcell*,field&)=0;
    virtual double ls1val(int,int,int)=0;
    virtual double ls2val(int,int,int)=0;
	virtual double ccipol_ls1val(lexer*,ghostcell*,double,double,double)=0;
	virtual double ccipol_ls2val(lexer*,ghostcell*,double,double,double)=0;
    virtual void ls1get(int,int,int,double)=0;
    virtual void ls2get(int,int,int,double)=0;

    virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&)=0;
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)=0;
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)=0;
};

#endif