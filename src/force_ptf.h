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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include<iostream>
#include<fstream>
#include"field4.h"
#include"fieldint4.h"

class lexer;
class fdm_ptf;
class ghostcell;

#ifndef FORCE_PTF_H_
#define FORCE_PTF_H_

using namespace std;

class force_ptf : public increment
{   
    
public:

    force_ptf(lexer*,fdm_ptf*,ghostcell*,int);
    virtual ~force_ptf();
    virtual void start(lexer*,fdm_ptf*,ghostcell*);
    virtual void ini(lexer*,fdm_ptf*,ghostcell*);
    
private:
    
	void print_force_ptf(lexer*,fdm_ptf*,ghostcell*);
    void print_ini_ptf(lexer*,fdm_ptf*,ghostcell*);
    
    // force ptf variabes
    double F_x_cell,F_y_cell,F_z_cell;
    double F_x_tot,F_y_tot,F_z_tot;
    
    // printing
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int force_ptfprintcount;
    ofstream fout;
    const int ID;
    
    // parallelisation
    double xstart,ystart,xend,yend,zstart,zend;
};

#endif