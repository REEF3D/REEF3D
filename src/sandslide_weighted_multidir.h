/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef SANDSLIDE_WEIGHTED_MULITDIR_H_
#define SANDSLIDE_WEIGHTED_MULITDIR_H_

#include"norm_vec.h"
#include"bedslope.h"
#include"slice4.h"
#include"sandslide.h"

class sediment_fdm;

using namespace std;

class sandslide_weighted_multidir :  public sandslide, public norm_vec, public bedslope
{
public:
    sandslide_weighted_multidir(lexer*);
    virtual ~sandslide_weighted_multidir();

	void start(lexer*,ghostcell*,sediment_fdm*) override;

private:
    void slide(lexer*,ghostcell*,sediment_fdm*);

    slice4 fh;
    
    int gcval_topo,count;

    double fac1, fac2;
    double dh,maxdh,maxdhs,dh_corr;
    double slide_dh,slide_dhs;
	double teta, alpha, beta, gamma;
    double phi;
    
    void compute_fh(lexer*, ghostcell*, sediment_fdm*);
    double compute_slope(lexer* p, slice&, int i, int j, int di, int dj);
    double compute_weight(double excess_slope, int method);
    
    double tan_phi;       // Tangent of angle of repose
    double relax;         // Relaxation factor (0 < relax <= 1)
    int maxiter;          // Maximum iterations
    double tol;           // Convergence tolerance
    int weight_method;    // Weighting scheme selector
}; 

#endif

