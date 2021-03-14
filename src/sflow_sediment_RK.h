/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;

using namespace std;

#ifndef SFLOW_SEDIMENT_RK_H_
#define SFLOW_SEDIMENT_RK_H_

class sflow_sediment_RK
{
public:

    virtual void step1(lexer*, fdm2D*, ghostcell*, slice&, slice&, double)=0;
    virtual void step2(lexer*, fdm2D*, ghostcell*, slice&, slice&, double)=0;
    virtual void step3(lexer*, fdm2D*, ghostcell*, slice&, slice&, double)=0;
    virtual void step4(lexer*, fdm2D*, ghostcell*, slice&, slice&, double)=0;
    
};

#endif
