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
Author: Hans Bihs
--------------------------------------------------------------------*/

#define HYPRE_COMPILATION
#ifdef  HYPRE_COMPILATION

#include"solver2D.h"
#include"increment.h"
#include"vec2D.h"

using namespace std;

#ifndef SFLOW_JACOBI_BLOCK_H_
#define SFLOW_JACOBI_BLOCK_H_

class sflow_jacobi_block : public solver2D, public increment
{
public:

	sflow_jacobi_block(lexer*,ghostcell*);
	virtual ~sflow_jacobi_block();
	virtual void start(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int, int, double, cpt2D&);
	virtual void solve(lexer*, ghostcell*, matrix2D&, vec2D&, vec2D&, int, int, int&, int, double, cpt2D&);
	virtual void setup(lexer*, ghostcell*,int, cpt2D&);
    
private:
    
    double res_calc(lexer*, ghostcell*, matrix2D&, vec2D&, vec2D&, cpt2D&);
    
    void fillxvec1(lexer*, slice&, vec2D&);
    void fillxvec2(lexer*, slice&, vec2D&);
    void fillxvec4(lexer*, slice&, vec2D&);
    
    void finalize(lexer*, slice&, vec2D&, int);
    

    int num_iterations;
    double final_res_norm,residual;
   
	int numiter,count,q,qn;
    const double epsi;
    
    int *sizeS,*range;


	int margin;

};

#endif

#endif

