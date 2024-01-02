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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

class mooring;
class net;
class ddweno_f_nug;
class lexer;
class fdm;
class fdm2D;
class ghostcell;
class field;
class reinidisc;
class mooring;
class net;
class vrans;
#include<vector>

using namespace std;

#ifndef SIXDOF_DF_BASE_H_
#define SIXDOF_DF_BASE_H_

class sixdof_df_base
{
public:

    virtual void start_forcing(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool)=0;
	
    virtual void start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&)=0;
	virtual void initialize(lexer*,fdm*,ghostcell*,vector<net*>&)=0;
    
    virtual void isource(lexer*,fdm*,ghostcell*)=0;
    virtual void jsource(lexer*,fdm*,ghostcell*)=0;
    virtual void ksource(lexer*,fdm*,ghostcell*)=0;
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*)=0;
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*)=0;

};

#endif
