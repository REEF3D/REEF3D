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

#include"sflow_turbulence.h"
#include"increment.h"
#include"slice4.h"
#include<fstream>

class fdm2D;
class lexer;
class ghostcell;

#ifndef SFLOW_TURB_IO_VOID_H_
#define SFLOW_TURB_IO_VOID_H_

using namespace std;

class sflow_turb_io_void : public sflow_turbulence, public increment
{

public:
    sflow_turb_io_void(lexer*);
	virtual ~sflow_turb_io_void();
    
    virtual void print_2D(lexer*, fdm2D*, ghostcell*,ofstream&);
    
    virtual void kinget(int,int,double);
    virtual void epsget(int,int,double);
    
    virtual double kinval(int,int);
    virtual double epsval(int,int);
    
	virtual void name_pvtp(lexer*, fdm2D*, ghostcell*,ofstream&);
    virtual void name_vtp(lexer*, fdm2D*, ghostcell*,ofstream&, int*, int &);
    
    virtual void offset_vtp(lexer*, fdm2D*, ghostcell*,ofstream&, int*, int &);
    

};

#endif
