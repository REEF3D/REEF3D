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
class sediment_fdm;

#include<fstream>
#include<sstream>
#include<vector>

using namespace std;



class sediment
{
public:

	virtual void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*)=0;
    virtual void ini_cfd(lexer*,fdm*,ghostcell*)=0;
    virtual void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*)=0;
    virtual void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*)=0;
    
    
    virtual void start_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*)=0;
    virtual void ini_nhflow(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void start_susp_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*)=0;
    virtual void update_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*)=0;
    
    virtual void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&)=0;
    virtual void ini_sflow(lexer*, fdm2D*, ghostcell*)=0;
    virtual void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*)=0;
    
    
    
    //
    void relax(lexer*,ghostcell*) override {};
	double bedshear_point(lexer*,ghostcell*) override {};
    
    double qbeval(int,int) override {};
    void qbeget(int,int,double) override {};
    
    double bedzhval(int,int) override {};

    void write_state_particles(lexer *, ofstream&) override {};
    void read_state_particles(lexer *, ifstream&) override {};
    
    void ctimesave(lexer*, fdm*) override {};
    
    virtual void print_probes(lexer*, ghostcell*,sediment_fdm*, ioflow*)=0;
    
    void print_2D_bedload(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_bedload(lexer*, ghostcell*,std::vector<char>&, size_t&) override {};
	void name_ParaView_parallel_bedload(lexer*,ofstream&) override {};
    void name_ParaView_bedload(lexer*, ostream&, int*, int &) override {};
    void offset_ParaView_2D_bedload(lexer*, int*, int &) override {};
    void offset_ParaView_bedload(lexer*, int*, int &) override {};
    
	void print_2D_bedshear(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_bedshear(lexer*, ghostcell*,std::vector<char>&, size_t&) override {};
	void name_ParaView_parallel_bedshear(lexer*,ofstream&) override {};
    void name_ParaView_bedshear(lexer*, ostream&, int*, int &) override {};
    void offset_ParaView_2D_bedshear(lexer*, int*, int &) override {};
    void offset_ParaView_bedshear(lexer*, int*, int &) override {};
    
    void print_2D_parameter1(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_parameter1(lexer*, ghostcell*,std::vector<char>&, size_t&) override {};
	void name_ParaView_parallel_parameter1(lexer*,ofstream&) override {};
    void name_ParaView_parameter1(lexer*, ostream&, int*, int &) override {};
    void offset_ParaView_2D_parameter1(lexer*, int*, int &) override {};
    void offset_ParaView_parameter1(lexer*, int*, int &) override {};
    
    void print_2D_parameter2(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_parameter2(lexer*, ghostcell*,std::vector<char>&, size_t&) override {};
	void name_ParaView_parallel_parameter2(lexer*,ofstream&) override {};
    void name_ParaView_parameter2(lexer*, ostream&, int*, int &) override {};
    void offset_ParaView_2D_parameter2(lexer*, int*, int &) override {};
    void offset_ParaView_parameter2(lexer*, int*, int &) override {};

};

#endif
