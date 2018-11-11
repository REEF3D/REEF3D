/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#define HYPRE_COMPILATION
#ifdef  HYPRE_COMPILATION

#include"solver2D.h"
#include"increment.h"
#include"vec2D.h"
#include"_hypre_utilities.h"
#include"HYPRE_struct_ls.h"

using namespace std;

#ifndef HYPRE_STRUCT2D_H_
#define HYPRE_STRUCT2D_H_

class hypre_struct2D : public solver2D, public increment
{
public:

	hypre_struct2D(lexer*,fdm2D*,ghostcell*);
	virtual ~hypre_struct2D();
	virtual void start(lexer*,fdm2D*, ghostcell*, slice&, vec2D&, vec2D&, int, int, double);
	virtual void solve(lexer*,fdm2D*, ghostcell*, vec2D&, vec2D&, int, int, int&, int, double, cpt2D&);
	virtual void setup(lexer*,fdm2D*, ghostcell*,int, cpt2D&);
    
    void fillxvec4(lexer*,fdm2D*,slice&);
	void fillbackvec(lexer*,fdm2D*,slice&,vec2D&,int);
    
    void make_grid(lexer*, ghostcell*);
	void fill_matrix(lexer*,fdm2D*, ghostcell*,slice&);
	
	void create_solvers(lexer*,ghostcell*);
    void delete_solvers(lexer*,ghostcell*);

private:
    
// HYPRE 
	HYPRE_StructGrid     grid;
	HYPRE_StructStencil  stencil;
	HYPRE_StructMatrix   A;
	HYPRE_StructVector   rhs;
	HYPRE_StructVector   x;
	HYPRE_StructSolver   solver;
	HYPRE_StructSolver   precond;

	int *ilower,*iupper;
    double *values;
    int num_iterations;
    double final_res_norm;
	int stencil_indices[7];
	int nentries;
   
	int numiter,count,q;

};

#endif

#endif

