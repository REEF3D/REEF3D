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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_motionext_file.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_motionext_file::read_format_1(lexer *p, ghostcell *pgc)
{
    char name[100];
	double val,val0,val1;
    double sign,beta,s;
	int count;
	
	sprintf(name,"6DOF_motion.dat");

    try {
        auto table = readTable(name, colnum);
        ptnum = table.size();
        
        // copy data directly
        data = std::move(table);
    } catch(const std::exception& e) {
        cout << "Error: " << e.what() << endl;
        pgc->final(EXIT_FAILURE);
    }
    
    ts = data[0][0];
    te = data[ptnum-1][0];
    
    //if(p->mpirank==0)
    //cout<<"6DOF_motion  ts: "<<ts<<" te: "<<te<<endl;
    
// add deltas
    for(qn=0;qn<ptnum;++qn)
    data[qn][0] += p->X241;
    
}
