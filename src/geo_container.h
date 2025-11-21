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

#ifndef GEO_CONTAINER_H_
#define GEO_CONTAINER_H_

#include"increment.h"

class lexer;
class ghostcell;

using namespace std;

class geo_container : increment
{
public:
	geo_container(lexer*);
	virtual ~geo_container();
    
    void create_obj(lexer*, int, int , int);
    
    void delete_obj(lexer*, int ID);

    
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
