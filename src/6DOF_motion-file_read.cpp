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

// open file and count
	ifstream file(name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no '6DOF_motion.dat' file found")<<endl<<endl;

    
    count=0;
	while(!file.eof())
	{
        for(qn=0;qn<colnum;++qn)
        file>>val;
    
	++count;
	}
	ptnum=count;
    
	file.close();
    
// allocate
    p->Darray(data,val,colnum);
    

// re.open file
    file.open (name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no '6DOF_motion.dat' file found")<<endl<<endl;
    
 // read file   
    rowcount=colcount=0;
	while(!file.eof())
	{
        for(qn=0;qn<colnum;++qn)
        file>>data[rowcount][qn];
        
        ++rowcount;
	}
    
    ts = data[0][0];
    te = data[ptnum-1][0];
    
// add deltas
    for(qn=0;qn<ptnum;++qn)
    data[qn][0] += p->X241;
    
}
