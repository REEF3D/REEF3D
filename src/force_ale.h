/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"nodefill.h"
#include"fieldint5.h"
#include"field5.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;

#ifndef FORCE_ALE_H_
#define FORCE_ALE_H_

using namespace std;

class force_ale :  public increment
{

public:
	force_ale(lexer*,fdm*,ghostcell*,int);
	virtual ~force_ale();
	virtual void start(lexer*,fdm*,ghostcell*);
    virtual void ini(lexer*,fdm*,ghostcell*);

private:	
	
    void force_ale_force(lexer*,fdm*,ghostcell*);
	void print_force_ale(lexer*,fdm*,ghostcell*);
    void print_ini(lexer*,fdm*,ghostcell*);
    
    
    // force_ale
    double Fx,Fy,Fz;	
    int is,js;
    const int ID;
    
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    
    int force_aleprintcount;
    
    ofstream fout;

};

#endif


