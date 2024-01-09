/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nsewave.h"


using namespace std;

#ifndef NSEWAVE_V_H_
#define NSEWAVE_V_H_

class nsewave_v : public nsewave
{
public:
    nsewave_v(lexer*, fdm*, ghostcell*,heat*&,concentration*&);
	virtual ~nsewave_v();
    
    virtual void start(lexer*, fdm*, ghostcell*, momentum*, diffusion*, turbulence*, convection*, 
                        pressure*, poisson*, solver*, solver*, ioflow*, vrans*, sixdof_df_base*, vector<net*>&);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*);
    void update(lexer*,fdm*,ghostcell*,slice&);

};

#endif
