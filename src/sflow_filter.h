/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"sflow_diffusion.h"
#include"increment.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;
class sliceint;

#ifndef SFLOW_EDIFF_H_
#define SLFOW_EDIFF_H_

using namespace std;

class sflow_filter : public increment
{
public:
	sflow_filter(lexer*);
	virtual ~sflow_filter();

    virtual void filter(lexer*, fdm2D*, ghostcell*);
	virtual void filter1(lexer*, fdm2D*, ghostcell*);
	virtual void filter2(lexer*, fdm2D*, ghostcell*);
    virtual void filter4(lexer*, fdm2D*, ghostcell*);
    
    slice1 f1x,f1y;
    slice2 f2x,f2y;
    slice4 f4x,f4y;

};

#endif
