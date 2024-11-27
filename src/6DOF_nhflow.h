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

#include"6DOF.h"
#include<vector>
#include"increment.h"
#include"6DOF_obj.h"

class lexer;
class fdm2D;
class fdm_nhf;
class ghostcell;
class net;
class slice;
class fdm;

using namespace std;

#ifndef SIXDOF_NHFLOW_H_
#define SIXDOF_NHFLOW_H_

class sixdof_nhflow : public sixdof, public increment
{
public:
	
    sixdof_nhflow(lexer*, ghostcell*);
    virtual ~sixdof_nhflow();
    
    virtual void start_cfd(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool);
    virtual void start_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*,vector<net*>&,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool);
    virtual void start_sflow(lexer*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool);
    
    void start_twoway(lexer*,fdm_nhf*,ghostcell*,vrans*,vector<net*>&,int,double*,double*,double*,slice&,slice&,bool);
    void start_oneway(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,slice&,slice&,bool);
    void start_shipwave(lexer*,fdm_nhf*,ghostcell*,int,bool);
       
	virtual void ini(lexer*,ghostcell*);
    virtual void initialize(lexer*, fdm*, ghostcell*, vector<net*>&);
    virtual void initialize(lexer*, fdm_nhf*, ghostcell*, vector<net*>&);
	
    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource(lexer*,fdm_nhf*,ghostcell*,slice&);
    virtual void jsource(lexer*,fdm_nhf*,ghostcell*,slice&);
    virtual void ksource(lexer*,fdm_nhf*,ghostcell*,slice&);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);
    
private:
	
    // hires gradient
    double limiter(double v1, double v2);
    
    double denom,val,r,phival;
    double dfdx_plus,dfdx_min,dfdy_plus,dfdy_min,dfdx,dfdy;


    int number6DOF;
    vector<sixdof_obj*> fb_obj;

    slice4 press;


};

#endif
