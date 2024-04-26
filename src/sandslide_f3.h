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

#include"norm_vec.h"
#include"bedslope.h"
#include"slice4.h"
#include"sandslide.h"

using namespace std;

#ifndef SANDSLIDE_F3_H_
#define SANDSLIDE_F3_H_

class sandslide_f3 :  public sandslide, public norm_vec, public bedslope
{
public:
    sandslide_f3(lexer*);
    virtual ~sandslide_f3();

	virtual void start(lexer*,ghostcell*, sediment_fdm*);

private:

    void slide(lexer*,ghostcell*, sediment_fdm*);
	
    slice4 fh;
    
    int gcval_topo,count;

    double fac1, fac2;
    double dh,maxdh,maxdhs,dh_corr;
    double slide_dh,slide_dhs;
	double teta, alpha, beta, gamma;
    double phi;
}; 

#endif

