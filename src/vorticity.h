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

class fdm;
class lexer;
class convection;
class diffusion;
class solver;
class ghostcell;
class ioflow;

#include<fstream>

#ifndef VORTICITY_H_
#define VORTICITY_H_

using namespace std;

class vorticity
{

public:

	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&)=0;


    virtual void name_pvtk(lexer*, fdm*, ghostcell*,ofstream&)=0;
    virtual void name_vtk(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)=0;
    virtual void offset_vtk(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)=0;
};

#endif


