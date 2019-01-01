/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_f::motion_vec(lexer *p, fdm* a, ghostcell* pgc)
{
    Uext=Vext=Wext=0.0;
    
    if(p->simtime>=ts_motion && p->simtime<=te_motion)
    {
	
    double vel = (motion[tcount_motion+1][1]-motion[tcount_motion][1])/(motion[tcount_motion+1][0]-motion[tcount_motion][0]);
	
    Uext = vel*nvecx;
	Vext = vel*nvecy;
	Wext = vel*nvecz;
    
	if(p->simtime>motion[tcount_motion+1][0])
	++tcount_motion;
    }
	
}

void sixdof_f::read_motionvec(lexer *p, fdm* a, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1;
    double sign,beta,s;
	int count,ptnum;
	
	sprintf(name,"motionfile_vec.dat");

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no 'motionfile_vec.dat' file found")<<endl<<endl;
	
	count=0;
	while(!file.eof())
	{
	file>>val0>>val1;
	++count;
	}
	
	file.close();

    
    ptnum=count;
	
	p->Darray(motion,ptnum,3);
	
	file.open ("motionfile_vec.dat", ios_base::in);
	
	count=0;
	while(!file.eof())
	{
        file>>val0>>val1;
        
        motion[count][0] = val0;
        motion[count][1] = val1;
        ++count;
	}
	
	ts_motion = motion[0][0];
	te_motion = motion[ptnum-1][0];
    
    
    // vec
    vecx = p->X221_xe - p->X221_xs;
    vecy = p->X221_ye - p->X221_ys;
    vecz = p->X221_ze - p->X221_zs;
    
    double length = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
    
    length = length>1.0e-20?length:1.0e20;
    
    nvecx = vecx/length;
    nvecy = vecy/length;
    nvecz = vecz/length;
    
    tcount_motion=0;
    
}