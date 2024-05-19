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
class fdm_fnpf;
class ghostcell;
class field;
class ioflow;
class wave_theory;
class slice;

using namespace std;

#ifndef FNPF_PRINT_WSFLINE_H_
#define FNPF_PRINT_WSFLINE_H_

class fnpf_print_wsfline : public boundarycheck
{
public:
    fnpf_print_wsfline(lexer*,fdm_fnpf*,ghostcell*);
	virtual ~fnpf_print_wsfline();

	void start(lexer*, fdm_fnpf*, ghostcell*,ioflow*,slice &f);


private:
    void ini_location(lexer*, fdm_fnpf*, ghostcell*);
    void sort(double*, double*, int*, int,int);
    void remove_multientry(lexer*,double*, double*, int*, int&);

    int conv(double);
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

