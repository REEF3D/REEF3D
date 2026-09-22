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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::read_format_piston(lexer *p, ghostcell *pgc)
{
    char name[100];
	double val,val0,val1;
    double sign,beta,s;
	int count,qn;
    
    colnum = 2;
	
    if(p->mpirank==0)
    cout<<"6DOF_motion piston wavemaker "<<endl;
    
	sprintf(name,"wavemaker.dat");

// open file and count
	ifstream file(name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no 'wavemaker.dat' file found")<<endl<<endl;
    
    count=0;
	while(!file.eof())
	{
	file>>val0>>val1;
	if(val0>=p->B117)
	++count;
	}
    ptnum=count;
	
	file.close();

// allocate
    p->Darray(kinematics,ptnum,colnum);
    

// re.open file
    file.open (name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no 'wavemaker.dat' file found")<<endl<<endl;
    
 // read file   
    rowcount=colcount=0;
	while(!file.eof())
	{
        for(qn=0;qn<colnum;++qn)
        file>>kinematics[rowcount][qn];
        
        //cout<<kinematics[rowcount][0]<<" "<<kinematics[rowcount][1]<<endl;
        
        ++rowcount;
	}
    
    ts = kinematics[0][0];
    te = kinematics[ptnum-1][0];
    
    // newline check
    while(te<1.0e-10 && ptnum>0)
    {
    --ptnum;
    te = kinematics[ptnum-1][0];
    }
    
    if(p->mpirank==0)
    cout<<"wavemaker.dat  ts: "<<ts<<" te: "<<te<<" ptnum: "<<ptnum<<endl;
    
// add deltas
    for(qn=0;qn<ptnum;++qn)
    kinematics[qn][0] += p->X241;
}

void sixdof_obj::read_format_flap(lexer *p, ghostcell* pgc)
{
}

void sixdof_obj::read_format_flap_double(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1,val2;
    double sign1,sign2;
    double beta,s,sign;
	int count;
	
	sprintf(name,"wavemaker.dat");

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	{
		cout<<endl<<("no 'wavemaker.dat' file found")<<endl<<endl;

	}
	
	count=0;
	while(!file.eof())
	{
	file>>val0>>val1>>val2;
	if(val0>=p->B117)
	++count;
	}
	
	file.close();

    
    ptnum=count;
	
	p->Darray(kinematics,ptnum,5);
	
	file.open ("wavemaker.dat", ios_base::in);
	
	count=0;
	while(!file.eof())
	{
	
        file>>val0>>val1>>val2;
        
        if(val0>=p->B117)
        {
        kinematics[count][0] = val0-p->B117;
        kinematics[count][1] = val1;
        kinematics[count][2] = val2;
        ++count;
        }
	}
	
	ts = kinematics[0][0];
	te = kinematics[ptnum-1][0];
    
    // convert from angles to X
    if(p->B116==2)
    for(int qn=0; qn<ptnum; ++qn)
    {
    sign1 = -(kinematics[qn][1]/(fabs(kinematics[qn][1])>1.0e-20?fabs(kinematics[qn][1]):1.0e20));
    sign2 = -(kinematics[qn][2]/(fabs(kinematics[qn][2])>1.0e-20?fabs(kinematics[qn][2]):1.0e20));
    
    kinematics[qn][1] = sign1*fabs(sin(kinematics[qn][1])*(p->X172_z2-p->X172_z1));
    kinematics[qn][2] = sign2*fabs(sin(kinematics[qn][2])*(p->X172_ze-p->X172_z2));
    }
    
    // calculate vertical component
    for(int qn=0; qn<ptnum; ++qn)
    {
    sign = -1.0;
    
    beta = asin(kinematics[qn][1]/(p->X172_z2-p->X172_z1));
    s = 2.0*(p->X172_z2-p->X172_z1)*sin(0.5*beta);
    kinematics[qn][3] = sign*sqrt(pow(s,2.0) - pow(kinematics[qn][1],2.0));
    
    beta = asin(kinematics[qn][2]/(p->X172_ze-p->X172_z2));
    s = 2.0*(p->X172_ze-p->X172_z2)*sin(0.5*beta);
    kinematics[qn][4] = sign*sqrt(pow(s,2.0) - pow(kinematics[qn][2],2.0));
    }
    
    /*
    cout<<"flap kinematics"<<endl;
    if(p->mpirank==0)
    for(int qn=0; qn<ptnum; ++qn)
    cout<<kinematics[qn][1]<<" "<<kinematics[qn][2]<<" | "<<kinematics[qn][3]<<" "<<kinematics[qn][4]<<endl;
    */
	    
}
