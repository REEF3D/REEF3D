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

void sflow_turb_io::print_2D(lexer *p, fdm2D *b, ghostcell *pgc, std::ostream &result)
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
    result<<"<PDataArray type=\"Float32\" Name=\"kin\"/>\n";
	
	if(p->A260==1)
	result<<"<PDataArray type=\"Float32\" Name=\"epsilon\"/>\n";
	if(p->A260==2)
    result<<"<PDataArray type=\"Float32\" Name=\"omega\"/>\n";
}

void sflow_turb_io::name_vtp(lexer *p, fdm2D *b, ghostcell *pgc,std::ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"kin\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	if(p->A260==1)
	result<<"<DataArray type=\"Float32\" Name=\"epsilon\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	if(p->A260==2)
    result<<"<DataArray type=\"Float32\" Name=\"omega\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}
    
void sflow_turb_io::offset_ParaView_2D(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sflow_turb_io::timesource(lexer* p, fdm2D *b, slice &fn)
{
    int count=0;
    SLICELOOP4
    {
        b->M.p[count] += 1.0/p->dt;

        b->rhsvec.V[count] += b->L(i,j) + fn(i,j)/p->dt;

	++count;
    }
}

void sflow_turb_io::wall_bc(lexer* p, fdm2D *b, slice &f)
{
    // SFLOW side walls are free-slip in the momentum equations (sflow_momentum_func): zero gradient of k and
    // epsilon / omega there and at dry neighbours, without wall functions
    n=0;
    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        if(p->flagslice4[Im1J]<0 || p->wet[Im1J]==0)
        {
        b->rhsvec.V[n] -= b->M.s[n]*f(i,j);
        b->M.s[n] = 0.0;
        }

        if(p->flagslice4[Ip1J]<0 || p->wet[Ip1J]==0)
        {
        b->rhsvec.V[n] -= b->M.n[n]*f(i,j);
        b->M.n[n] = 0.0;
        }

        if(p->flagslice4[IJm1]<0 || p->wet[IJm1]==0)
        {
        b->rhsvec.V[n] -= b->M.e[n]*f(i,j);
        b->M.e[n] = 0.0;
        }

        if(p->flagslice4[IJp1]<0 || p->wet[IJp1]==0)
        {
        b->rhsvec.V[n] -= b->M.w[n]*f(i,j);
        b->M.w[n] = 0.0;
        }
        }

    ++n;
    }

    // dry cells: identity rows, f = 0
    n=0;
    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
        b->M.p[n] = 1.0;
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.w[n] = 0.0;
        b->M.e[n] = 0.0;
        b->rhsvec.V[n] = 0.0;
        }
    ++n;
    }
}
