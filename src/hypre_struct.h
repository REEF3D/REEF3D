/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"solver.h"
#include"increment.h"
#include"vec.h"
#include"_hypre_utilities.h"
#include"HYPRE_sstruct_ls.h"
 
using namespace std;

#ifndef HYPRE_STRUCT_H_
#define HYPRE_STRUCT_H_

class hypre_struct : public solver, public increment
{
public:

	hypre_struct(lexer*,ghostcell*,int,int);
	virtual ~hypre_struct();
    
	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, int);
    virtual void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int);
    virtual void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    virtual void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    virtual void startM(lexer*, ghostcell*, double*, double*, double*, int);
    
    void start_solver1234(lexer*,fdm*, ghostcell*, field&, vec&,int);
    void start_solver4f(lexer*, ghostcell*, field&, vec&, matrix_diag&, int);
    void start_solver4V(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void start_solver5(lexer*,fdm*, ghostcell*, field&, vec&, int);
    void start_solver7(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void start_solver8(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void start_solver9(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    
    void solve(lexer*,ghostcell*);
    void solve1234(lexer*);
    
	void fillxvec1(lexer*,fdm*,field&);
    void fillxvec2(lexer*,fdm*,field&);
    void fillxvec3(lexer*,fdm*,field&);
    void fillxvec4(lexer*,fdm*,field&);
    
    
    
    void make_grid(lexer*, ghostcell*);
    void make_grid_2Dvert(lexer*, ghostcell*);
    void make_grid_15pt(lexer*, ghostcell*);
    void make_grid_2D_9pt(lexer*, ghostcell*);
    
    void fill_matrix1(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix1_2Dvert(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix2(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix2_2Dvert(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix3(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix3_2Dvert(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix4(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix4_2Dvert(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix4V(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix4V_2D(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix4f(lexer*, ghostcell*,field&, vec&, matrix_diag&);
    void fill_matrix4f_2Dvert(lexer*, ghostcell*,field&, vec&, matrix_diag&);
    void fill_matrix7(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix7_2Dvert(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix8(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix8_2Dvert(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix9(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix9_2Dvert(lexer*, ghostcell*,double*, vec&, matrix_diag&);

    virtual void fillbackvec1(lexer*,field&,int);
    virtual void fillbackvec2(lexer*,field&,int);
    virtual void fillbackvec3(lexer*,field&,int);
    virtual void fillbackvec4(lexer*,field&,int);
    virtual void fillbackvec4V(lexer*,double*,int);
    
    virtual void fillbackvec7(lexer*,double*,int);
    virtual void fillbackvec8(lexer*,double*,int);
    virtual void fillbackvec9(lexer*,double*,int);
	
	void create_solver1234(lexer*,ghostcell*);
    void delete_solver1234(lexer*,ghostcell*);

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
    double *values;
    int num_iterations;
    double final_res_norm;
	int stencil_indices[7];
    int periodic[3];
	int nentries;
   
	int numiter,count,q;
    
    const int solve_type,precon_type;
    
    
    int *CVAL4;

};

#endif

#endif

