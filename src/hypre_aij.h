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
#include"HYPRE_krylov.h"
#include"HYPRE.h"
#include"HYPRE_parcsr_ls.h"

using namespace std;

#ifndef HYPRE_AIJ_H_
#define HYPRE_AIJ_H_

class hypre_aij : public solver, public increment
{
public:

	hypre_aij(lexer*,fdm*,ghostcell*);
	virtual ~hypre_aij();
	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, int);
    virtual void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int);
    virtual void startM(lexer*, ghostcell*, double*, double*, double*, int);
    virtual void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    
	void solve(lexer*,fdm*, ghostcell*, vec&, vec&, int, int, int&);

	void fillxvec1(lexer*,fdm*,field&);
    void fillxvec2(lexer*,fdm*,field&);
    void fillxvec3(lexer*,fdm*,field&);
    void fillxvec4(lexer*,fdm*,field&);
	void fillbackvec(lexer*,fdm*,field&,vec&,int);
    
    
    void make_grid(lexer*,ghostcell*);
    void delete_grid(lexer*,ghostcell*);
	void fill_matrix_7p(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix_13p(lexer*,fdm*, ghostcell*,field&);
    void fill_matrix_19p(lexer*,fdm*, ghostcell*,field&);
    
    void fillbackvec_F(lexer*,double*,double*,int);
    void fillbackvec_F_v2(lexer*,double*,double*,int);
    
    
    void create_solvers(lexer*,ghostcell*);
    void delete_solvers(lexer*,ghostcell*);
    
    // FNPF Laplace solver 
    virtual void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void make_grid_F(lexer*, ghostcell*);
    void fill_matrix_F_7p(lexer*,ghostcell*, matrix_diag&,double*,double*,vec&);
    void fill_matrix_F_13p(lexer*, ghostcell*, matrix_diag&,double*,double*,vec&);
    void fill_matrix_F_19p(lexer*, ghostcell*, matrix_diag&,double*,double*,vec&);
    void fill_matrix_F_7p_v2(lexer*,ghostcell*, matrix_diag&,double*,double*,vec&);
    void fill_matrix_F_13p_v2(lexer*, ghostcell*, matrix_diag&,double*,double*,vec&);
    void fill_matrix_F_19p_v2(lexer*, ghostcell*, matrix_diag&,double*,double*,vec&);
	
    
    
private:
    
// HYPRE 
    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;
    
    HYPRE_Solver solver, precond;
    
    vec xvec;
   
	double val[19];
	int col[19];
	int *rows;
	int rownum;
    int num_iterations;
    double final_res_norm;
   
    int is,ie,js,je,ks,ke;
	
    
// -------------
	
	int *sizeM;

	int numiter,count,q;
	double resi,y,residual;
	double p1,p2,p3;
	int margin;
	int matlength;

};

#endif

#endif

