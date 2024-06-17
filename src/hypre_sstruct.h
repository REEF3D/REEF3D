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


#ifndef HYPRE_SSTRUCT_H_
#define HYPRE_SSTRUCT_H_

#define HYPRE_COMPILATION

#ifdef HYPRE_COMPILATION

#include"solver.h"
#include"increment.h"
#include"vec.h"
#include"_hypre_utilities.h"
#include"HYPRE_sstruct_ls.h"
#include"HYPRE_parcsr_ls.h"
#include"HYPRE_krylov.h"
#include"HYPRE.h"

using namespace std;

class hypre_sstruct : public solver, public increment
{
public:

	hypre_sstruct(lexer*,fdm*,ghostcell*);
	virtual ~hypre_sstruct();
    
	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, int);
    virtual void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int);
    virtual void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    virtual void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    virtual void startM(lexer*, ghostcell*, double*, double*, double*, int);
    
    void start_solver1234(lexer*,fdm*, ghostcell*, field&, vec&, int);
    void start_solver5(lexer*,fdm*, ghostcell*, field&, vec&,int);
    void start_solver7(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void start_solver8(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void start_solver10(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void start_solverM(lexer*, ghostcell*, double*, double*, double*);
    
    void solve(lexer*);
    void solve1234(lexer*);
    
	void fillxvec1(lexer*,fdm*,field&);
    void fillxvec2(lexer*,fdm*,field&);
    void fillxvec3(lexer*,fdm*,field&);
    void fillxvec4(lexer*,fdm*,field&);
    
    void make_grid_7p(lexer*,fdm*, ghostcell*);
    void make_grid_13p(lexer*,fdm*, ghostcell*);
    void make_grid_15p(lexer*,fdm*, ghostcell*);
    void make_grid_2Dvert_9p(lexer*,fdm*, ghostcell*);
    
    void fill_matrix1(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix2(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix3(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix4(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix7(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix8(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    void fill_matrix10(lexer*, ghostcell*,double*, vec&, matrix_diag&);
    
    void fill_matrixM(lexer*, ghostcell*,double*, double*, double*);
    void fill_matrixM_2Dvert(lexer*, ghostcell*,double*, double*, double*);

    void fillbackvec1(lexer*,field&,int);
    void fillbackvec2(lexer*,field&,int);
    void fillbackvec3(lexer*,field&,int);
    void fillbackvec4(lexer*,field&,int);

    void fillbackvec7(lexer*,double*,int);
    void fillbackvec8(lexer*,double*,int);
    void fillbackvec10(lexer*,double*,int);
    
    void fillbackvecM(lexer*,double*);
    
    
	
	void create_solver1234(lexer*,ghostcell*);
    void delete_solver1234(lexer*,ghostcell*);

    void create_solver5(lexer*,ghostcell*);
    void delete_solver5(lexer*,ghostcell*);
    

private:
    
//  HYPRE 
    HYPRE_SStructGrid     grid;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructStencil  stencil;
    HYPRE_SStructMatrix   A;
    HYPRE_SStructVector   b;
    HYPRE_SStructVector   x;
    HYPRE_SStructSolver   solver;
    HYPRE_SStructSolver   precond;
    HYPRE_Solver solver_csr, precond_csr;
    HYPRE_SStructVariable vartypes[1];
    
    HYPRE_ParCSRMatrix    par_A;
    HYPRE_ParVector       par_b;
    HYPRE_ParVector       par_x;
   
    int kend;
    int numparts;
    int part;
    int dimensions;
    int variable;
    int numvar;
    int object_type;
   

	int *ilower,*iupper;
    double *values;
    int num_iterations;
    double final_res_norm;
	int stencil_indices[13];
	int nentries;
   
	int numiter,count,q;

};

#endif

#endif

