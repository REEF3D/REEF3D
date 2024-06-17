/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"printer_nhflow.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void printer_nhflow::parallelData(lexer *p, ghostcell* pgc)
{	
	int num=0;
    
    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	
    outputFormat->parallelFileName(name, "NHFLOW", num);


	ofstream result;
	result.open(name);
    if(result.is_open())
	{
        outputFormat->beginningParallel(p,result);
        
        if(p->P16==1)
        {
        result<<"<FieldData>"<<endl;
        result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
        result<<"</DataArray>"<<endl;
        result<<"</FieldData>"<<endl;
        }
        
        result<<"<PPointData>"<<endl;
        result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
        result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;
        result<<"<PDataArray type=\"Float32\" Name=\"omega_sig\"/>"<<endl;
        result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
        if(p->P23==1)
        result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
        if(p->P25==1)
        result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;
        result<<"</PPointData>"<<endl;
        
        outputFormat->endingParallel(result,"NHFLOW",p->M10,num);

        result.close();
    }
    else
    cout<<"Failed to open output file."<<endl;
}