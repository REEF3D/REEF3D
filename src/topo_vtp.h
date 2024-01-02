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

#include"increment.h"
#include"fieldint5.h"
#include"field5.h"

class lexer;
class fdm;
class ghostcell;
class sediment;

#ifndef TOPO_VTP_H_
#define TOPO_VTP_H_

using namespace std;

class topo_vtp :  public increment
{

public:
	topo_vtp(lexer*,fdm*,ghostcell*);
	virtual ~topo_vtp();
	virtual void start(lexer*,fdm*,ghostcell*,sediment*);

private:
	void print(lexer*,fdm*,ghostcell*,sediment*);
    void pvtp(lexer*,fdm*,ghostcell*,sediment*);
    void header(lexer*,fdm*,ghostcell*);
    
    void name_iter(lexer*,fdm*,ghostcell*);
    void name_time(lexer*,fdm*,ghostcell*);
    void piecename(lexer*,fdm*,ghostcell*,int);

    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    double ddn;
    int gcval_phi;
    double printtime,printtime2;
	int topoprintcount;
	int polygon_sum,polygon_num;

};

#endif


