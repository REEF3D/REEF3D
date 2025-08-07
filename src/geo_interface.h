/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef GEO_INTERFACE_H_
#define GEO_INTERFACE_H_

#include"increment.h"

class lexer;
class ghostcell;
class fdm;
class fdm_nhf;
class fdm2D;

using namespace std;

class geo_interface : increment
{
public:
	geo_interface(lexer*);
	virtual ~geo_interface();
    
    void create_geometry(lexer*);
    
    void delete_geometry(lexer*);
    
    // CFD
    void solid_update_cfd(lexer*,fdm*);
    void fb_update_cfd(lexer*,fdm*);
    void topo_update_cfd(lexer*,fdm*);
    void por_update_cfd(lexer*,fdm*);

    
private:
    double **tri;
    int trinum;
    int ID;
    int type;
    double xs,xe,ys,ye,zs,ze;
   
   
    int q,iin;
    float ffn;
    int offset[100];
};

#endif
