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

#ifndef SANDSLIDE_STEEPEST_DESCENT_H_
#define SANDSLIDE_STEEPEST_DESCENT_H_

#include"norm_vec.h"
#include"bedslope.h"
#include"slice4.h"
#include"sandslide.h"

class sediment_fdm;

using namespace std;

class sandslide_steepest_descent :  public sandslide, public norm_vec, public bedslope
{
public:
    sandslide_steepest_descent(lexer*);
    virtual ~sandslide_steepest_descent();

	void start(lexer*,ghostcell*,sediment_fdm*) override final;

private:
    void find_steepest_neighbor(lexer* p, slice& topo, int i, int j,
                                 int& i_steep, int& j_steep, 
                                 double& max_slope, double& dist_steep);

    void slide(lexer*,ghostcell*,sediment_fdm*);

    slice4 fh;
    
    int gcval_topo,count;

    double fac1, fac2;
    double dh,maxdh,maxdhs,dh_corr;
    double slide_dh,slide_dhs;
	double teta, alpha, beta, gamma;
    double phi;
    
    double relax;
}; 

#endif

