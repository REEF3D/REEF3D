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

#include"sflow_turb_io.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_turb_io::sflow_turb_io(lexer* p) : kin(p), eps(p)
{

}

sflow_turb_io::~sflow_turb_io()
{
}

void sflow_turb_io::print_2D(lexer *p, fdm2D *b, ghostcell *pgc, ofstream &result)
{
    iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(kin));
	result.write((char*)&ffn, sizeof (float));
	}
    
    
    iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(eps));
	result.write((char*)&ffn, sizeof (float));
	}
    
}

void sflow_turb_io::kinget(int ii, int jj, double val)
{
    kin(ii,jj)=val;
}

void sflow_turb_io::epsget(int ii, int jj, double val)
{
    eps(ii,jj)=val;
}
    
double sflow_turb_io::kinval(int ii, int jj)
{
    val=kin(ii,jj);

    return val;
}
    
double sflow_turb_io::epsval(int ii, int jj)
{
    val=eps(ii,jj);

    return val;
}
    
void sflow_turb_io::name_pvtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"kin\"/>"<<endl;
	
	if(p->A260==1)
	result<<"<PDataArray type=\"Float32\" Name=\"epsilon\"/>"<<endl;
	if(p->A260==2 || p->A260==5)
    result<<"<PDataArray type=\"Float32\" Name=\"omega\"/>"<<endl;
}

void sflow_turb_io::name_vtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"kin\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	if(p->A260==1)
	result<<"<DataArray type=\"Float32\" Name=\"epsilon\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	if(p->A260==2 || p->A260==5)
    result<<"<DataArray type=\"Float32\" Name=\"omega\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}
    
void sflow_turb_io::offset_vtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}