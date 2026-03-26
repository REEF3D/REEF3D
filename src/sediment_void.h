/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SEDIMENT_VOID_H_
#define SEDIMENT_VOID_H_

#include"sediment.h"

using namespace std;

class sediment_void : public sediment
{
public:
    sediment_void();
	virtual ~sediment_void();
    
    void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*) override final;
    void ini_cfd(lexer*,fdm*,ghostcell*) override final;
    void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*) override final;
    void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*) override final;
    
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*) override final;
    void ini_nhflow(lexer*, fdm_nhf*, ghostcell*) override final;
    void start_susp_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*) override final;
    void update_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override final;
    
    void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&) override final;
    void ini_sflow(lexer*, fdm2D*, ghostcell*) override final;
    void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*) override final;
    
    // ---
	
    void relax(lexer*,ghostcell*) override final;
	double bedshear_point(lexer*,ghostcell*) override final;
    
    double qbeval(int,int) override final;
    void qbeget(int,int,double) override final;
    
    double bedzhval(int,int) override final;
    
    void ctimesave(lexer*, fdm*) override final;
    
    void print_probes(lexer*, ghostcell*,sediment_fdm*, ioflow*) override final {};
    void print_particles(lexer*,sediment_fdm*) override final {};
    
    void print_2D_bedload(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_bedload(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_bedload(lexer*,ofstream&) override final;
    void name_ParaView_bedload(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_bedload(lexer*, int*, int &) override final;
    void offset_ParaView_bedload(lexer*, int*, int &) override final;
    
	void print_2D_bedshear(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_bedshear(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_bedshear(lexer*,ofstream&) override final;
    void name_ParaView_bedshear(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_bedshear(lexer*,int*, int &) override final;
    void offset_ParaView_bedshear(lexer*, int*, int &) override final;
    
    void print_2D_parameter1(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_parameter1(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_parameter1(lexer*,ofstream&) override final;
    void name_ParaView_parameter1(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_parameter1(lexer*, int*, int &) override final;
    void offset_ParaView_parameter1(lexer*, int*, int &) override final;
    
    void print_2D_parameter2(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_parameter2(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_parameter2(lexer*,ofstream&) override final;
    void name_ParaView_parameter2(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_parameter2(lexer*, int*, int &) override final;
    void offset_ParaView_parameter2(lexer*, int*, int &) override final;
    
    void print_3D_CPM(lexer*, ghostcell*,  std::vector<char>&, size_t&) override final {};
    void name_ParaView_parallel_CPM(lexer*, ofstream&) override final {};
    void name_ParaView_CPM(lexer*, ostream&, int*, int &) override final {};
    void offset_ParaView_CPM(lexer*, int*, int &) override final {};
};

#endif
