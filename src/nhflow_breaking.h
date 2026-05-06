/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef NHFLOW_BREAKING_H_
#define NHFLOW_BREAKING_H_

#include"nhflow_fsf.h"
#include"sliceint4.h"
#include"slice4.h"

using namespace std;

class nhflow_breaking : public increment 
{
public:
	nhflow_breaking(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_breaking();
    
    void breaking(lexer*,fdm_nhf*,ghostcell*,slice&,slice&,double);
    
    void breaking_baquet(lexer*,fdm_nhf*,ghostcell*,slice&,slice&,double);

    void filter(lexer*, fdm_nhf*, ghostcell*, slice&);

    
private:

    double visc;
    int count_n;
    
    sliceint4 bx,by,brkflag;
};

#endif
