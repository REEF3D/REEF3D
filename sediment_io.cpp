/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bedshear.h"

double sediment_f::bedshear_point(lexer *p, fdm *a,ghostcell *pgc)
{
	double tau_eff,shearvel_eff,shields_eff;
	
	pbedshear->taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);

	return tau_eff;
}

void sediment_f::fill_bedk(lexer *p, fdm *a,ghostcell *pgc)
{
    SLICELOOP4
    a->bedk(i,j)=0;
    
    SLICELOOP4
    KLOOP
    PBASECHECK
    if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
    a->bedk(i,j)=k;
    
    
    SLICELOOP1
    {
    k=a->bedk(i,j);
    a->P(i,j) = a->u(i,j,k);
    }
    
    SLICELOOP2
    {
    k=a->bedk(i,j);
    a->Q(i,j) = a->v(i,j,k);
    }
    
    pgc->gcsl_start1(p,a->P,10);
	pgc->gcsl_start2(p,a->Q,11);
}

void sediment_f::fill_bss(lexer *p, fdm *a,ghostcell *pgc)
{
	double tau_eff,shearvel_eff,shields_eff;

    SLICELOOP4
    {
	pbedshear->taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
	bedtau(i,j) = tau_eff;
	}
	
	pgc->gcsl_start4(p,bedtau,1);
    pgc->dgcslpol(p,bedtau,p->dgcsl4,p->dgcsl4_count,14);
    bedtau.ggcpol(p);
}

void sediment_f::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
    
    ALOOP
    bss(i,j,k) = bedtau(i,j);
	
    pgc->start4(p,bss,1);
    pgc->dgcpol(p,bss,p->dgc4,p->dgc4_count,14);
    bss.ggcpol(p);
	
	iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
	ffn=float(p->ipol4_a(bss));
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
	ffn=float(p->ccipol4_a(bss,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));
	}
}

void sediment_f::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"bedshear\"/>"<<endl;
}

void sediment_f::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"bedshear\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
}