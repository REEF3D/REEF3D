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

#ifndef FNPF_PRINT_WSF_H_
#define FNPF_PRINT_WSF_H_

#include"increment.h"
#include"amr_grid_lookup.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice;

using namespace std;

class fnpf_print_wsf : public increment
{
public:
    // the mesh refinement: the points read the finest grid (amr_grid_lookup.h), nullptr: level 0
    const amr_column_lookup *amr = nullptr;

    fnpf_print_wsf(lexer*,fdm_fnpf*);
	virtual ~fnpf_print_wsf();

	void height_gauge(lexer*, fdm_fnpf*, ghostcell*,slice&);


private:
    void ini_location(lexer*, fdm_fnpf*);
    void fill_eta(lexer*, fdm_fnpf*, ghostcell*,slice&);
    void fill_deta(lexer*, fdm_fnpf*, ghostcell*,slice&);
    void fill_Uhorz(lexer*, fdm_fnpf*, ghostcell*,slice&);
    
	double *x,*y;
	int gauge_num;

    int *iloc,*jloc,*flag;
    double *wsf,*deta,*Uhorz;
    int n;
    ofstream wsfout,detaout,Uhorzout;
public:
    void flush_files() { wsfout.flush(); detaout.flush(); Uhorzout.flush(); }   // a rebalance stop (G 46)
private:
    const int fileFlushMaxCount;
    

};

#endif
