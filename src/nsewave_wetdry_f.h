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

#include"nsewave_wetdry.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef NSEWAVE_WETDRY_F_H_
#define NSEWAVE_WETDRY_F_H_

class nsewave_wetdry_f : public nsewave_wetdry
{
public:
    nsewave_wetdry_f(lexer*, fdm*, ghostcell*);
	virtual ~nsewave_wetdry_f();
    
    virtual void start(lexer*, fdm*, ghostcell*);
    virtual void ini(lexer*, fdm*, ghostcell*);

};

#endif