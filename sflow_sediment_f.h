/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_sediment.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;

using namespace std;

#ifndef SFLOW_SEDIMENT_F_H_
#define SFLOW_SEDIMENT_F_H_

class sflow_sediment_f : public sflow_sediment, public increment
{
public:

    sflow_sediment_f(lexer*, fdm2D*);
	virtual ~sflow_sediment_f();

	virtual void ini(lexer*, fdm2D*, ghostcell*);
    virtual void start(lexer*, fdm2D*, ghostcell*, slice&, slice&,slice&);
    
private:
    void sediment_algorithm(lexer*, fdm2D*, ghostcell*,slice&,slice&,slice&);
    
    void bedslope(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void bedshear(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void shields(lexer*, fdm2D*, ghostcell*);
    void bedshear_slope(lexer*, fdm2D*, ghostcell*);
    
    void bedload(lexer*, fdm2D*, ghostcell*);
    void bedload_vanRijn(lexer*, fdm2D*, ghostcell*);
    
    void exner(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&);
    void sandslide(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    
    
    slice4 tau,taucr,alpha,teta,gamma,phi;
    
    double starttime;
    
};

#endif