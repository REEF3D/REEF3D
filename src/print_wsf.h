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

#ifndef PRINT_WSF_H_
#define PRINT_WSF_H_

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

class amr_grid_lookup;

class print_wsf : public boundarycheck
{
public:
    void flush_files() { wsfout.flush(); }   // a rebalance stop (G 46)
    print_wsf(lexer*,fdm*,ghostcell*,int);
	virtual ~print_wsf();

	void height_gauge(lexer*, fdm*, ghostcell*,field&);

    const amr_grid_lookup *amr = nullptr;     // CFD mesh refinement: the free surface on the finest grid


private:
    void ini_location(lexer*, fdm*, ghostcell*);
    void amr_gauge(lexer*, int);
    void write(lexer*, fdm*, ghostcell*);
	
	double *x,*y;
	int gauge_num;

    int *iloc,*jloc,*flag;
    double *wsf;
    int n;
    ofstream wsfout;
    const int fileFlushMaxCount;

};

#endif
