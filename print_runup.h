/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
class ioflow;
class wave_theory;

using namespace std;

#ifndef PRINT_RUNUP_H_
#define PRINT_RUNUP_H_

class print_runup  : public boundarycheck
{
public:
    print_runup(lexer*,fdm*,ghostcell*);
	virtual ~print_runup();

	void start(lexer*, fdm*, ghostcell*);
	void runup(lexer*, fdm*, ghostcell*);
	void cone_ini(lexer*, fdm*, ghostcell*);
	void cone(lexer*, fdm*, ghostcell*);
	void print(lexer*, fdm*, ghostcell*);
	
	double **line;
	int line_num,line_num_all;
	
	double *runup_x,*runup_y,*runup_z;
	int *runup_active;
	
	
	int cut_count;
	double *cut_x,*cut_y,*cut_z;
	int *cut_ID,*cut_active;
	
	int *cutall_count,cuttotal_count;
	double *cutall_x,*cutall_y,*cutall_z;
	int *cutall_ID;
	
	int *displ;
	
	
	ofstream maxout,fsfout;
	


private:
	char name[100];
	int num,n;

};

#endif

