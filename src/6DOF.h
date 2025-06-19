/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef SIXDOF_H_
#define SIXDOF_H_

#include<vector>

class lexer;
class fdm;
class fdm_nhf;
class fdm2D;
class ghostcell;
class vrans;
class net;
class field;
class slice;

using namespace std;

class sixdof
{
public:
    virtual void start_cfd(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool)=0;
    virtual void start_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*,vector<net*>&,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool)=0;
    virtual void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool)=0;
    
    virtual void ini(lexer*,ghostcell*)=0;
    virtual void initialize(lexer*, fdm*, ghostcell*, vector<net*>&)=0;
    virtual void initialize(lexer*, fdm2D*, ghostcell*, vector<net*>&)=0;
    virtual void initialize(lexer*, fdm_nhf*, ghostcell*, vector<net*>&)=0;
	
    virtual void isource(lexer*,fdm*,ghostcell*)=0;
    virtual void jsource(lexer*,fdm*,ghostcell*)=0;
    virtual void ksource(lexer*,fdm*,ghostcell*)=0;
    
    virtual void isource(lexer*,fdm_nhf*,ghostcell*,slice&)=0;
    virtual void jsource(lexer*,fdm_nhf*,ghostcell*,slice&)=0;
    virtual void ksource(lexer*,fdm_nhf*,ghostcell*,slice&)=0;
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*)=0;
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*)=0;
};

#endif
