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
Authora: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef PARTRES2_H_
#define PARTRES2_H_

#include"increment.h"
#include"slice4.h"

class lexer;
class fdm;
class ghostcell;
class sediment_fdm;
class turbulence;
class part;

using namespace std;


class partres2 : public increment
{
public:
    partres2(lexer *);
    ~partres2();
    
    void move_RK2(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);
    
    void advec_plain(lexer*, fdm*, part&, sediment_fdm*, turbulence*, 
                        double*, double*, double*, double*, double*, double*, 
                        double&, double&, double&, double);
    void advec_pic(lexer*, fdm*, part&, sediment_fdm*, turbulence*, 
                        double*, double*, double*, double*, double*, double*, 
                        double&, double&, double&, double);
    
    // drag
    double drag_model(lexer *, double, double, double);
    double drag_coefficient(double);
    
    void timestep(lexer*, ghostcell*, part*);
    
    
    // relax
    void relax_ini(lexer*);
    void relax(lexer*, ghostcell*, sediment_fdm*);
    double rf(lexer*, double, double);     
    double r1(lexer*, double, double);
    double distcalc(lexer*, double , double, double , double, double);
        
        
private:
    const int irand;
	const double drand;
    
    
    
};

#endif