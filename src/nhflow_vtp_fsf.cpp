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

#include"nhflow_vtp_fsf.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_vtp_fsf::nhflow_vtp_fsf(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	if(p->I40==0)
    {
	p->printtime=0.0;
    }
	
	printcount=0;
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_NHFLOW_VTP_FSF",0777);
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
}

nhflow_vtp_fsf::~nhflow_vtp_fsf()
{
}

void nhflow_vtp_fsf::start(lexer *p, fdm_nhf *d, ghostcell* pgc)
{	
    print2D(p,d,pgc);
}

void nhflow_vtp_fsf::print2D(lexer *p, fdm_nhf *d, ghostcell* pgc)
{	
    //pgc->gcsl_start4(p,d->eta,gcval_eta);
    //pgc->gcsl_start4(p,d->Fifsf,gcval_fifsf);
    
    pgc->gcsl_start4(p,d->test2D,1);
    
    SLICELOOP4
    {
    if(d->breaking(i,j)>=1)
    d->breaking_print(i,j)=double(d->breaking(i,j));
        
    if(d->breaking(i,j)==0)
    d->breaking_print(i,j)=0.0;   
    }
    
    //pgd->gcsl_start4(p,d->breaking_print,50);
    
    d->eta.ggcpol(p);
    
	if(p->mpirank==0)
    pvtu(p,d,pgc);
    
	name_iter(p,d,pgc);
	
	
	// Open File
	ofstream result;
	result.open(name, ios::binary);
    
    // offsets
    n=0;
	offset[n]=0;
	++n;
	
	// Points
    offset[n]=offset[n-1]+8*(p->pointnum2D)*3+4;
    ++n;
	
	// velocity
	offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
	++n;
    
    // Fifsf
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	
	// WL
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // breaking
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // coastline
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // wetdry
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // test
    if(p->P23==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }
    
    // Hs
    if(p->P110==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }
	
	// Cells
    offset[n]=offset[n-1] + 4*p->polygon_sum*3+4;
    ++n;
    offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;
	offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;
	
	
	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum2D<<"\" NumberOfPolys=\""<<p->polygon_sum<<"\">"<<endl;
    
    n=0;
	result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float64\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
	
    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"eta\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"detadt\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"WL\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"coastline\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"wetdry\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    if(p->P23==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    if(p->P110==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"Hs\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    result<<"</PointData>"<<endl;

    

    result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
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

    //ddn=float(p->sl_ipol4(d->eta) + p->wd);
    
    //ddn=float(0.5*(d->eta(i,j) + d->eta(i,j))  + p->wd);
    ddn=p->sl_ipol4eta(p->wet,d->eta, d->bed)+p->wd;
	result.write((char*)&ddn, sizeof (double));
	}
	
    //  Velocities
    iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
    k = p->knoz-1;
    
	if(p->j_dir==0)
    {
	jj=j;
    j=0;
	ffn=float(d->U[IJK]);
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.5*(d->U[IJK]+d->U[IJp1K]));
    
	result.write((char*)&ffn, sizeof (float));


	if(p->j_dir==0)
    {
	jj=j;
    j=0;
	ffn=float(d->UH[IJK]);
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.5*(d->V[IJK]+d->V[IJp1K]));
    
	result.write((char*)&ffn, sizeof (float));


	if(p->j_dir==0)
    {
	jj=j;
    j=0;
	ffn=float(d->W[IJK]);
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.5*(d->W[IJK]+d->W[IJp1K]));
    
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Eta
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4eta_wd(p->wet,d->eta,d->bed));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Detadt
    k=p->knoz-1;
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(d->detadt));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  WL
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->WL));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Breaking
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
        
	ffn=float(p->sl_ipol4(d->breaking_print));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Coastline
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->coastline));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Wetdry
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
    ffn = 0.25*float((p->wet[IJ]+p->wet[Ip1J]+p->wet[IJp1]+p->wet[Ip1Jp1]));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  test
    if(p->P23==1)
    {
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->test2D));
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    //  Hs
    if(p->P110==1)
    {
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->Hs));
	result.write((char*)&ffn, sizeof (float));
	}
    }

    //  Connectivity
    iin=4*(p->polygon_sum)*3;
    result.write((char*)&iin, sizeof (int));
    SLICEBASELOOP
	{
	// Triangle 1
	iin=int(d->nodeval2D(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->nodeval2D(i,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->nodeval2D(i,j))-1;
	result.write((char*)&iin, sizeof (int));
	
	
	// Triangle 2
	iin=int(d->nodeval2D(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->nodeval2D(i,j))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->nodeval2D(i-1,j))-1;
	result.write((char*)&iin, sizeof (int));
	}
    
    
    //  Offset of Connectivity
    iin=4*(p->polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=(n+1)*3;
	result.write((char*)&iin, sizeof (int));
	}
    
//  Cell types
    iin=4*(p->polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

    result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
	
	++printcount;

}


