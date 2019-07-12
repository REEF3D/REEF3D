/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fnpf_sg.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

#ifndef GRID_SIGMA_H_
#define GRID_SIGMA_H_

class grid_sigma : public increment
{
public:
	grid_sigma(lexer*);
	virtual ~grid_sigma();
    
    virtual void sigma_coord_ini(lexer*);
    virtual void sigma_ini(lexer*, fdm*, ghostcell*, slice&);
    virtual void sigma_update(lexer*, fdm*, ghostcell*, slice&);
    
    double sigmax(lexer*,field&,int);
    double sigmay(lexer*,field&,int);
    double sigmaz(lexer*,field&,int);

    
    slice4 Ex,Ey,Bx,By;
    slice4 Exx,Eyy,Bxx,Byy;
        
private:
    double sig;

};

#endif
