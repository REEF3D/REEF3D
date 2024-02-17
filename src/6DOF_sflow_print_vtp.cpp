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

#include<sys/stat.h>
#include<iostream>
#include<fstream>
#include"6DOF_sflow.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_sflow::print_ini_vtp(lexer *p, ghostcell *pgc)
{
	if(p->mpirank==0)
    {
        mkdir("./REEF3D_SFLOW_6DOF_VTP", 0777);
        mkdir("./REEF3D_SFLOW_6DOF", 0777);
    }
	
    ofstream print;
    char str[1000];

    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
	
    print.open(str);
	print<<"time \t XG \t YG \t ZG \t Phi \t Theta \t Psi"<<endl;
	print.close();
    
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
	
    print.open(str);
	print<<"time \t Ue [m/s] \t Ve [m/s] \t We [m/s] \t Pe [rad/s] \t Qe [rad/s] \t Re [rad/s]"<<endl;
    print.close();
}


void sixdof_sflow::print_vtp(lexer *p, ghostcell *pgc)
{
	int num=0;

	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;

    
    if(p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  || (p->simtime>printtime && p->P30>0.0)   || p->count==0))
    {
        printtime+=p->P30;
        
 
        char path[300];
        
       sprintf(path,"./REEF3D_SFLOW_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp",n6DOF,num);

        ofstream result;
        result.open(path, ios::binary);
  
        
    // ---------------------------------------------------
    n=0;

	offset[n]=0;
	++n;

    offset[n]=offset[n-1]+4*tricount*3*3 + 4;
    ++n;
    offset[n]=offset[n-1]+4*tricount*3 + 4;
    ++n;
    offset[n]=offset[n-1]+4*tricount*3 + 4;
    ++n;
	//---------------------------------------------

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<tricount*3<<"\" NumberOfPolys=\""<<tricount<<"\">"<<endl;

    n=0;
    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;

    result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;

	result<<"</Polys>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</PolyData>"<<endl;

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";


//  XYZ
	iin=4*tricount*3*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<tricount;++n)
	for(q=0;q<3;++q)
	{
	ffn=tri_x[n][q];
	result.write((char*)&ffn, sizeof (float));

	ffn=tri_y[n][q];
	result.write((char*)&ffn, sizeof (float));

	ffn=tri_z[n][q];
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Connectivity POLYGON
	int count=0;
    iin=4*tricount*3;
    result.write((char*)&iin, sizeof (int));
    for(n=0;n<tricount;++n)
	for(q=0;q<3;++q)
	{
	iin=count;
	result.write((char*)&iin, sizeof (int));
	++count;
	}

//  Offset of Connectivity
    iin=4*tricount;
    result.write((char*)&iin, sizeof (int));
	iin=0;
	for(n=0;n<tricount;++n)
	{
	iin+= 3;//a->polygon_offset[n];
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*tricount;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<tricount;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();	

        ++p->printcount_sixdof;	
    }
}



