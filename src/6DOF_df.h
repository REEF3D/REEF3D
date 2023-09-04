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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_df_base.h"
#include"6DOF_df_object.h"
#include<vector>

class mooring;
class net;
class ddweno_f_nug;

using namespace std;

#ifndef SIXDOF_DF_H_
#define SIXDOF_DF_H_

class sixdof_df : public sixdof_df_base, public increment
{
public:
	sixdof_df(lexer*, fdm*, ghostcell*);
	virtual ~sixdof_df();
	virtual void start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&);
	virtual void initialize(lexer*,fdm*,ghostcell*,vector<net*>&);
    
    void start_forcing(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool);
    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);

private:
   
    int number6DOF;
    vector<sixdof_df_object*> p_df_obj;
    double alpha[3],gamma[3],zeta[3];

};

#endif
