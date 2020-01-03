/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"norm_vec.h"
#include"bedslope.h"
#include"field4a.h"
#include"sandslide.h"

using namespace std;

#ifndef SANDSLIDE_PDE_H_
#define SANDSLIDE_PDE_H_

class sandslide_pde :  public sandslide, public norm_vec, public bedslope
{
public:
    sandslide_pde(lexer*);
    virtual ~sandslide_pde();

	virtual void start(lexer*, fdm*,ghostcell*);

private:

    void slide(lexer*, fdm*,ghostcell*);
    void diff_update(lexer*, fdm*,ghostcell*);
    void topo_zh_update(lexer*,fdm*,ghostcell*);
	
    field4a fh,ci;
    
    int gcval_topo,count;

    double fac1, fac2;
    double dh,maxdh,maxdhs,dxs,dh_corr;
    double slide_dh,slide_dhs;
	double teta, alpha, beta, gamma;
    double phi;
}; 

#endif

