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

#include"concentration.h"
#include"increment.h"
#include"field4.h"
#include"fluid_update.h"
#include<fstream>

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef CONCENTRATION_IO_H_
#define CONCENTRATION_IO_H_

class concentration_io : public concentration, increment
{
public:
    concentration_io(lexer*,fdm*);
	virtual ~concentration_io();

    virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void ini(lexer*, fdm*, ghostcell*, concentration *pconcentration);
    virtual double val(int,int,int);

    virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);

    field4 C;
	fluid_update *pupdate;

private:

	float ffn;
	double ddn;
	int n,iin;
	
	double fx(double,double,double,double,double);
	double fz(double,double,double,double,double);
};

#endif

