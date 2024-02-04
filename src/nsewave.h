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

#include<vector>

class convection;
class pressure;
class solver;
class fdm;
class lexer;
class ghostcell;
class fluid_update;
class heat;
class concentration;
class ioflow;
class slice;
class momentum;
class diffusion;
class poisson;
class vrans;
class turbulence;
class net;
class sixdof;

using namespace std;

#ifndef NSEWAVE_H_
#define NSEWAVE_H_

class nsewave 
{
public:    
    virtual void start(lexer*, fdm*, ghostcell*, momentum*, diffusion*, turbulence*, convection*, pressure*, 
                       poisson*, solver*, solver*, ioflow*, vrans*, sixdof*, vector<net*>&)=0;
    
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*)=0;

        

};

#endif
