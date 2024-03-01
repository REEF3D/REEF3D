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
Author: Hans Bihs
--------------------------------------------------------------------*/

#define HYPRE_COMPILATION

#ifdef HYPRE_COMPILATION

#include"solver_fnpf.h"
#include"increment.h"
#include"vec.h"
#include"fieldint4.h"
#include"_hypre_utilities.h"
#include"HYPRE_sstruct_ls.h"

using namespace std;

#ifndef HYPRE_STRUCT_FNPF_H_
#define HYPRE_STRUCT_FNPF_H_

class hypre_struct_fnpf : public solver_fnpf, public increment
{
public:

	hypre_struct_fnpf(lexer*,ghostcell*,int,int);
	virtual ~hypre_struct_fnpf();

    virtual void startF(lexer*, ghostcell*, double*, double*, double*, int);
    
    void start_solver8(lexer*, ghostcell*, double*, double*, double*);
    
    virtual void solve(lexer*,ghostcell*);
    
    void make_grid(lexer*, ghostcell*);
    void make_grid_2Dvert(lexer*, ghostcell*);

    void fill_matrix8(lexer*, ghostcell*, double*, double*, double*);
    void fill_matrix8_2Dvert(lexer*, ghostcell*, double*, double*, double*);


    virtual void fillbackvec8(lexer*,double*,double*,double*);
	

    void create_solver5(lexer*,ghostcell*);
    void delete_solver5(lexer*,ghostcell*);
    

private:
    
// HYPRE 
   HYPRE_StructGrid     grid;
   HYPRE_StructStencil  stencil;
   HYPRE_SStructGraph   graph;
   HYPRE_StructMatrix   A;
   HYPRE_StructVector   b;
   HYPRE_StructVector   x;
   HYPRE_StructSolver   solver;
   HYPRE_StructSolver   precond;
   

	int *ilower,*iupper;
    int num_iterations;
    double final_res_norm;
	int stencil_indices[15];
	int nentries;
   
	int numiter,count,q;
    
    const int solve_type,precon_type;


};

#endif

#endif

