/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef SFLOW_BICGSTAB_H_
#define SFLOW_BICGSTAB_H_

#include "solver2D.h"
#include "increment.h"

#include <vector>

using namespace std;

class sflow_bicgstab final : public solver2D, public increment
{
public:
    sflow_bicgstab(lexer*,ghostcell*);
    virtual ~sflow_bicgstab() = default;
    void start(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int) override final;

private:
    void solve(lexer*, ghostcell*, matrix2D&, vec2D&, vec2D&, int, int&);

    void fillxvec(lexer*,slice&,vec2D&);
    void finalize(lexer*,slice&);

    double res_calc(lexer*, matrix2D&, ghostcell*, std::vector<double>&);
    void matvec_axb(lexer*, matrix2D&, std::vector<double>&, std::vector<double>&);
    void matvec_std(lexer*, matrix2D&, std::vector<double>&, std::vector<double>&);

    void precon_setup(lexer*, matrix2D&,ghostcell*);
    void precon_solve(lexer*,ghostcell*,std::vector<double>&,std::vector<double>&);

    std::vector<double> sj,rj,r0,vj,tj,pj,ph,sh,aii,x,rhs;

    int *flagslice;

    double final_res_norm;
    double alpha,beta,w1,w2,w,residual,norm_vj,norm_r0,norm_sj,norm_rj;
    double r_j1,r_j,sigma;

    int count,q;

    int ulast,vlast,wlast;

    static constexpr double epsi = 1.0e-19;
};

#endif
