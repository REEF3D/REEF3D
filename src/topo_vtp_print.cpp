/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
	psed->offset_ParaView_2D_bedload(p,offset,n);

    // sediment parameters 1
	if(p->P77==1)
	psed->offset_ParaView_2D_parameter1(p,offset,n);

    // sediment parameters 2
	if(p->P78==1)
	psed->offset_ParaView_2D_parameter2(p,offset,n);

    // bed shear stress
	if(p->P79>=1)
	psed->offset_ParaView_2D_bedshear(p,offset,n);
    
	//End Data
    
    offset[n]=offset[n-1] + 4*polygon_sum*3 + 4;
    ++n;
    offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
	offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
	//---------------------------------------------
	
	

	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	result<<"<PolyData>\n";
	result<<"<Piece NumberOfPoints=\""<<p->pointnum2D<<"\" NumberOfPolys=\""<<polygon_sum<<"\">\n";

    n=0;
    result<<"<Points>\n";
    result<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";
	
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    
    if(p->P76==1)
	psed->name_ParaView_bedload(p,pgc,result,offset,n);
    
    if(p->P77==1)
	psed->name_ParaView_parameter1(p,pgc,result,offset,n);

    if(p->P78==1)
	psed->name_ParaView_parameter2(p,pgc,result,offset,n);

	if(p->P79>=1)
	psed->name_ParaView_bedshear(p,pgc,result,offset,n);
    result<<"</PointData>\n";

    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	++n;
    result<<"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	result<<"</Polys>\n";
	

    result<<"</Piece>\n";
    result<<"</PolyData>\n";

//----------------------------------------------------------------------------

    result<<"<AppendedData encoding=\"raw\">\n_";

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

	result<<"\n</AppendedData>\n";
    result<<"</VTKFile>\n";

	result.close();	
}







