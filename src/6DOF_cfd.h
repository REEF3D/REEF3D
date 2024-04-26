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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF.h"
#include"6DOF_obj.h"
#include<vector>

class mooring;
class net;
class ddweno_f_nug;

using namespace std;

#ifndef SIXDOF_CFD_H_
#define SIXDOF_CFD_H_

class sixdof_cfd : public sixdof, public increment
{
public:
	sixdof_cfd(lexer*, fdm*, ghostcell*);
	virtual ~sixdof_cfd();

    void start_twoway(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool);
    void start_twoway(lexer*,fdm_nhf*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool);
    virtual void start_oneway(lexer*,ghostcell*,slice&);
    
    virtual void ini(lexer*,ghostcell*);
    virtual void initialize(lexer*, fdm*, ghostcell*, vector<net*>&);
    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource(lexer*,fdm_nhf*,ghostcell*,slice&);
    virtual void jsource(lexer*,fdm_nhf*,ghostcell*,slice&);
    virtual void ksource(lexer*,fdm_nhf*,ghostcell*,slice&);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);

private:
   void setup(lexer*,fdm*,ghostcell*);
   
    int number6DOF;
    vector<sixdof_obj*> fb_obj;

};

#endif
