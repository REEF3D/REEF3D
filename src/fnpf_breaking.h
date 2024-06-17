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

#ifndef FNPF_BREAKING_H_
#define FNPF_BREAKING_H_

#include"fnpf_fsf.h"
#include"sliceint4.h"

class fnpf_laplace;
class field;
class fnpf_convection;
class fnpf_ddx;
class fnpf_etadisc;
class fnpf_coastline;
class solver2D;

using namespace std;

class fnpf_breaking : public fnpf_fsf, public increment 
{
public:
	fnpf_breaking(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_breaking();
    
    virtual void breaking_algorithm(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&,slice&,double);

    
    
    void filter(lexer*, fdm_fnpf*,ghostcell*, slice&);

    double ivel,jvel,kvel;
    
private:
    double rb3(lexer*,double);
    double rb4(lexer*,double);
    
    double dist3,dist4,db;
    
    double visc;
    
    sliceint4 bx,by;
    int count_n;
    
};

#endif
