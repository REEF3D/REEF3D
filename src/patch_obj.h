/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

#ifndef PATCH_OBJ_H_
#define PATCH_OBJ_H_

class patch_obj : public increment
{
public:
	patch_obj(lexer*);
	virtual ~patch_obj();
    
    void patch_obj_ini(lexer *p, ghostcell *pgc);
    
    void patch_obj_gcb_generate(lexer *p, ghostcell *pgc);
    
    // Patch DATA
    int ID;
    int IO;
    int gcb_count;
    int **gcb;
    
    int Q_flag;
    double Q;
    
    int velocity_flag;
    double velocity;
    
    int pressure_flag;
    double pressure;
    
    int waterlevel_flag;
    double waterlevel;
    
        
private:
    
    

};

#endif
