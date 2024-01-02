/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef BEDPROBE_LINE_X_H_
#define BEDPROBE_LINE_X_H_

class bedprobe_line_x : public boundarycheck
{
public:
    bedprobe_line_x(lexer*,fdm*,ghostcell*);
	virtual ~bedprobe_line_x();

	void start(lexer*, fdm*, ghostcell*,ioflow*);


private:
    void ini_location(lexer*, fdm*, ghostcell*);
    void sort(double*, double*, int*, int,int);
    void remove_multientry(lexer*,double*, double*, int*, int&);

    int *jloc,**flag,**flag_all,*rowflag,*wsfpoints;
    double **wsf,**wsf_all;
    double **xloc, **xloc_all;
    double *yloc;
    int n,q;
    ofstream wsfout;

    double xcoor;
	
	wave_theory *pwave;

    int maxknox,sumknox;

};

#endif

