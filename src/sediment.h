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

#ifndef SEDIMENT_H_
#define SEDIMENT_H_

class lexer;
class fdm;
class fdm2D;
class fdm_nhf;
class convection;
class ghostcell;
class ioflow;
class reinitopo;
class suspended;
class topo;
class reinitopo;
class field;
class slice;
class solver;

#include<fstream>

using namespace std;



class sediment
{
public:

	virtual void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*)=0;
    virtual void ini_cfd(lexer*,fdm*,ghostcell*)=0;
    virtual void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*)=0;
    virtual void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*)=0;
    
    virtual void start_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*)=0;
    virtual void ini_nhflow(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void update_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*)=0;
    
    virtual void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&)=0;
    virtual void ini_sflow(lexer*, fdm2D*, ghostcell*)=0;
    virtual void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*)=0;
    
    
    
    //
    virtual void relax(lexer*,ghostcell*){};
	virtual double bedshear_point(lexer*,ghostcell*){};
    
    virtual double qbeval(int,int){};
    virtual void qbeget(int,int,double){};
    
    virtual double bedzhval(int,int){};

    virtual void write_state_particles(lexer *, ofstream&){};
    virtual void read_state_particles(lexer *, ifstream&){};
    
    virtual void ctimesave(lexer*, fdm*){};
    
    virtual void print_2D_bedload(lexer*, ghostcell*,ofstream&){};
    virtual void print_3D_bedload(lexer*, ghostcell*,ofstream&){};
	virtual void name_pvtu_bedload(lexer*, ghostcell*,ofstream&){};
    virtual void name_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtp_bedload(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &){};
    
	virtual void print_2D_bedshear(lexer*, ghostcell*,ofstream&){};
    virtual void print_3D_bedshear(lexer*, ghostcell*,ofstream&){};
	virtual void name_pvtu_bedshear(lexer*, ghostcell*,ofstream&){};
    virtual void name_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtp_bedshear(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &){};
    
    virtual void print_2D_parameter1(lexer*, ghostcell*,ofstream&){};
    virtual void print_3D_parameter1(lexer*, ghostcell*,ofstream&){};
	virtual void name_pvtu_parameter1(lexer*, ghostcell*,ofstream&){};
    virtual void name_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtp_parameter1(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &){};
    
    virtual void print_2D_parameter2(lexer*, ghostcell*,ofstream&){};
    virtual void print_3D_parameter2(lexer*, ghostcell*,ofstream&){};
	virtual void name_pvtu_parameter2(lexer*, ghostcell*,ofstream&){};
    virtual void name_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtp_parameter2(lexer*, ghostcell*,ofstream&, int*, int &){};
    virtual void offset_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &){};

};

#endif
