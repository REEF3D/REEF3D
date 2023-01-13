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

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;
class turbulence;

using namespace std;

#ifndef FLOWFILE_OUT_H_
#define FLOWFILE_OUT_H_

class flowfile_out : public boundarycheck
{
public:
    flowfile_out(lexer*,fdm*,ghostcell*);
	virtual ~flowfile_out();

	void start(lexer*, fdm*, ghostcell*,turbulence*);


private:
    void filename(lexer*,fdm*,ghostcell*);
    void header_file(lexer*, fdm*, ghostcell*);
    void header_file_ini(lexer*, fdm*, ghostcell*);
    
    void initialize(lexer*, fdm*, ghostcell*);
    void ini_location(lexer*, fdm*, ghostcell*);
    
    void write_data(lexer*, fdm*, ghostcell*);

    
    
	char name[450];
    char headername[450];

    int **flag;
    double **U,**V,**W,**P,**K,**E,**VT,**LS;
    int n,q;
    int count;
    int elnum;
    int *iloc;
	const int probenum;
	const double eps;
    ofstream *fileout;
    ofstream *headerout;

	double xloc,yloc,zloc;
	double xp,yp,zp;

	int filecount;
    
    double ddn;
    float ffn;
    int iin;
    int Ni,Nj,Nk;
    
};

#endif
