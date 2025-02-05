/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"increment.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class slice;
class sixdof;
class vrans;
class mooring;
class net;
class fsi;
class nhflow_reinidisc_fsf;
#include<vector>

using namespace std;

#ifndef SFLOW_FORCING_H_
#define SFLOW_FORCING_H_

class sflow_forcing : public increment
{
public:
	sflow_forcing(lexer*);
	virtual ~sflow_forcing();
    
    void forcing(lexer*, fdm2D*, ghostcell*, sixdof *p6dof, 
                 int, double, slice&, slice&, slice&, slice&, slice&, bool);
    
    void forcing_ini(lexer*, fdm2D*, ghostcell*);
    
private:
    slice1 fx;
    slice2 fy;
    slice4 fz;
    
    double uf, vf, wf;
    int forcing_flag;

 
};

#endif