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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;

#ifndef print_porous_H_
#define print_porous_H_

using namespace std;

class print_porous :  public increment
{

public:
	print_porous(lexer*,fdm*,ghostcell*);
	virtual ~print_porous();
	virtual void start(lexer*,fdm*,ghostcell*);
    virtual void print_vtp(lexer*,fdm*,ghostcell*);
	virtual void objects(lexer*,fdm*,ghostcell*);
	
    void box(lexer*,fdm*,ghostcell*,int);
    void cylinder_z(lexer*,fdm*,ghostcell*,int);
	void wedge_x(lexer*,fdm*,ghostcell*,int);
    void wedge_y(lexer*,fdm*,ghostcell*,int);
    void plate_x(lexer*,fdm*,ghostcell*,int);
    
    void box_veg(lexer*,fdm*,ghostcell*,int);
    void wedge_x_veg(lexer*,fdm*,ghostcell*,int);
    void wedge_y_veg(lexer*,fdm*,ghostcell*,int);

private:
    char name[100];
    int offset[200];
	int **polygon,*numvert;
    double **vertice;
	
    float ffn;
    int polygon_num,vertice_num,polygon_sum ,iin,q;
	double xs,xe,ys,ye,zs,ze;
	int vertice_alloc, polygon_alloc;


};

#endif


