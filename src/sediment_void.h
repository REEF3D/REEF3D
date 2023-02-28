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

#include"sediment.h"

using namespace std;

#ifndef SEDIMENT_VOID_H_
#define SEDIMENT_VOID_H_

class sediment_void : public sediment
{
public:
    sediment_void();
	virtual ~sediment_void();
    
    virtual void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*);
    virtual void ini_cfd(lexer*,fdm*,ghostcell*);
    virtual void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*);
    virtual void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*);
    
    virtual void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&);
    virtual void ini_sflow(lexer*, fdm2D*, ghostcell*);
    virtual void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*);
    
    // ---
	
    virtual void relax(lexer*,ghostcell*);
	virtual double bedshear_point(lexer*,fdm*,ghostcell*);
    
    virtual double qbeval(int,int);
    virtual void qbeget(int,int,double);
    
    virtual double bedzhval(int,int);
    
    virtual void ctimesave(lexer*, fdm*);
    
    virtual void print_2D_bedload(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_bedload(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedload(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_bedload(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &);
    
	virtual void print_2D_bedshear(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_bedshear(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedshear(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_2D_parameter1(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_parameter1(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameter1(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_parameter1(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_2D_parameter2(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_parameter2(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameter2(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_parameter2(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &);
};

#endif
