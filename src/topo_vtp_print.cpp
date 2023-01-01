/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"topo_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"

void topo_vtp::print(lexer* p, fdm* a, ghostcell *pgc, sediment *psed)
{
	if(p->mpirank==0)
    pvtp(p,a,pgc,psed);
	
	name_iter(p,a,pgc);

	ofstream result;
	result.open(name, ios::binary);

	//---------------------------------------------
    n=0;
	offset[n]=0;
	++n;
    // Points
    offset[n]=offset[n-1] + 8*p->pointnum2D*3 + 4;
    ++n;
    
	//Velocity
	offset[n]=offset[n-1] + 4*p->pointnum2D*3+ 4;
    ++n;
    // Elevation
	offset[n]=offset[n-1] + 4*p->pointnum2D + 4;
    ++n;
    
    // sediment bedlaod
	if(p->P76==1)
	psed->offset_vtp_bedload(p,pgc,result,offset,n);

    // sediment parameters 1
	if(p->P77==1)
	psed->offset_vtp_parameter1(p,pgc,result,offset,n);

    // sediment parameters 2
	if(p->P78==1)
	psed->offset_vtp_parameter2(p,pgc,result,offset,n);

    // bed shear stress
	if(p->P79>=1)
	psed->offset_vtp_bedshear(p,pgc,result,offset,n);
    
	//End Data
    
    offset[n]=offset[n-1] + 4*polygon_sum*3 + 4;
    ++n;
    offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
	offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
	//---------------------------------------------
	
	

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum2D<<"\" NumberOfPolys=\""<<polygon_sum<<"\">"<<endl;

    n=0;
    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float64\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
    if(p->P76==1)
	psed->name_vtu_bedload(p,pgc,result,offset,n);
    
    if(p->P77==1)
	psed->name_vtu_parameter1(p,pgc,result,offset,n);

    if(p->P78==1)
	psed->name_vtu_parameter2(p,pgc,result,offset,n);

	if(p->P79>=1)
	psed->name_vtu_bedshear(p,pgc,result,offset,n);
    result<<"</PointData>"<<endl;

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
	iin=8*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
    ddn=p->XN[IP1];
	result.write((char*)&ddn, sizeof (double));

	ddn=p->YN[JP1];
	result.write((char*)&ddn, sizeof (double));

	ddn=p->sl_ipol4(a->bed);
    result.write((char*)&ddn, sizeof (double));
	}
    
//  Velocities
    iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol1a(a->P));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->sl_ipol2a(a->Q));
	result.write((char*)&ffn, sizeof (float));

	ffn=0.0;
	result.write((char*)&ffn, sizeof (float));
	}
	
//  Elevation
	iin=4*p->pointnum2D;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(a->bed));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  sediment bedload
	if(p->P76==1)
    psed->print_2D_bedload(p,pgc,result);
    
    //  sediment parameter 1
	if(p->P77==1)
    psed->print_2D_parameter1(p,pgc,result);

    //  sediment parameter 2
	if(p->P78==1)
    psed->print_2D_parameter2(p,pgc,result);

    //  bed shear stress
	if(p->P79>=1)
    psed->print_2D_bedshear(p,pgc,result);
    

//  Connectivity
    iin=4*(polygon_sum)*3;
    result.write((char*)&iin, sizeof (int));
    SLICEBASELOOP
	{
	// Triangle 1
	iin=int(a->nodeval2D(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval2D(i,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval2D(i,j))-1;
	result.write((char*)&iin, sizeof (int));


	// Triangle 2
	iin=int(a->nodeval2D(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval2D(i,j))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval2D(i-1,j))-1;
	result.write((char*)&iin, sizeof (int));
	}


    // Offset of Connectivity
    iin=4*(polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<polygon_sum;++n)
	{
	iin=(n+1)*3;
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*(polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<polygon_sum;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();	
}







