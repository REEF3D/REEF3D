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

#ifndef SFLOW_PRINT_PROBE_DA_H_
#define SFLOW_PRINT_PROBE_DA_H_

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm2D;
class ghostcell;
class field;
class turbulence;

using namespace std;

class sflow_print_probe_da : public boundarycheck
{
public:
    sflow_print_probe_da(lexer*,fdm2D*,ghostcell*);
	virtual ~sflow_print_probe_da();

	void start(lexer*, fdm2D*, ghostcell*);


private:
    void ini_location(lexer*, fdm2D*, ghostcell*);
    void write(lexer*, fdm2D*, ghostcell*);
	char name[100];

    int *iloc,*jloc,*flag;
    int n,q;
	const int probenum;
    ofstream *pout;
	
	double uval,vval,wval,pval,eval;

};

#endif
