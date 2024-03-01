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

#include"multiphase_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reini.h"
#include"print_wsf.h"
	
void multiphase_f::print_3D(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
	pgc->dgcpol(p,ls1,p->dgc4,p->dgc4_count,14);
    ls1.ggcpol(p);
    pgc->dgcpol(p,ls2,p->dgc4,p->dgc4_count,14);
    ls2.ggcpol(p);
	pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
    a->ro.ggcpol(p);

    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
	ffn=float(p->ipol4(ls1));
	result.write((char*)&ffn, sizeof (float));
	}


	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

	TPLOOP
	{
	ffn=float(p->ipol4(ls2));
	result.write((char*)&ffn, sizeof (float));
	}

	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

	TPLOOP
	{
	ffn=float(p->ipol4_a(a->ro));
	result.write((char*)&ffn, sizeof (float));
	}
	
}

void multiphase_f::print_file(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->P351>0)
	pwsf1->height_gauge(p,a,pgc,ls1);
	
	if(p->P352>0)
    pwsf2->height_gauge(p,a,pgc,ls2);	
}

void multiphase_f::nodefill(lexer *p, fdm *a, ghostcell *pgc, field &eta)
{
	
}

double multiphase_f::ls1val(int ii, int jj, int kk)
{
	double val;
    i=ii;
    j=jj;
    k=kk;

    val=ls1(i,j,k);

    return val;
	
}

double multiphase_f::ls2val(int ii, int jj, int kk)
{
	double val;
    i=ii;
    j=jj;
    k=kk;

    val=ls2(i,j,k);

    return val;
	
}

double multiphase_f::ccipol_ls1val(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
	double val;

    val=p->ccipol4( ls1, xp, yp, zp);

    return val;
}

double multiphase_f::ccipol_ls2val(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
	double val;

    val=p->ccipol4( ls2, xp, yp, zp);

    return val;
}

void multiphase_f::ls1get(int ii, int jj, int kk, double val)
{
	i=ii;
    j=jj;
    k=kk;

    ls1(i,j,k)=val;
}

void multiphase_f::ls2get(int ii, int jj, int kk, double val)
{
	i=ii;
    j=jj;
    k=kk;

    ls2(i,j,k)=val;
}

void multiphase_f::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
	result<<"<PDataArray type=\"Float32\" Name=\"ls1\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ls2\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"rho\"/>"<<endl;
}

void multiphase_f::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
	result<<"<DataArray type=\"Float32\" Name=\"ls1\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"ls2\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"rho\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void multiphase_f::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}


