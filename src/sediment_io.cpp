/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"bedshear.h"

double sediment_f::bedshear_point(lexer *p, fdm *a,ghostcell *pgc)
{
	double tau_eff = s->tau_eff(i,j);
    
	return tau_eff;
}

void sediment_f::fill_bedk(lexer *p, fdm *a,ghostcell *pgc)
{
    double zval,xip,yip;
    SLICELOOP4
    a->bedk(i,j)=0;
    
    SLICELOOP4
    KLOOP
    PBASECHECK
    if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
    a->bedk(i,j)=k+1;
    
    
    SLICELOOP1
    {
    k=a->bedk(i,j);
    
    xip= p->XN[IP1];
	yip= p->YP[JP];
    zval = 0.5*(a->bedzh(i,j)+a->bedzh(i+1,j)) + 1.6*p->DZN[k];
    
    a->P(i,j) = p->ccipol1(a->u,xip,yip,zval);
    }
    
    SLICELOOP2
    {
    k=a->bedk(i,j);
    
    xip= p->XP[IP];
	yip= p->YN[JP1];
    zval = 0.5*(a->bedzh(i,j)+a->bedzh(i,j+1)) + 1.6*p->DZN[k];
    
    a->Q(i,j) = p->ccipol2(a->v,xip,yip,zval);
    }
    
    pgc->gcsl_start1(p,a->P,10);
	pgc->gcsl_start2(p,a->Q,11);
}

void sediment_f::print_3D_bedshear(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
    double tau_eff,shearvel_eff,shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
    
    if(p->P79==1)
    {
        
	// tau_eff
    pgc->gcsl_start4(p,s->tau_eff,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->tau_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // tau_crit
    pgc->gcsl_start4(p,s->tau_crit,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->tau_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    
    }
    
    
    if(p->P79==2)
    {
        
	// shearvel_eff
    pgc->gcsl_start4(p,s->shearvel_eff,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // shearvel_crit
    pgc->gcsl_start4(p,s->shearvel_crit,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    
    }
    
    
    if(p->P79==3)
    {
        
	// shields_eff
    pgc->gcsl_start4(p,s->shields_eff,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shields_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // shields_crit
    pgc->gcsl_start4(p,s->shields_crit,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shields_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    
    }
}

void sediment_f::name_pvtu_bedshear(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    if(p->P79==1)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"tau_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"tau_crit\"/>"<<endl;
    }
    
    if(p->P79==2)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"shearvel_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"shearvel_crit\"/>"<<endl;
    }
    
    if(p->P79==2)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"shields_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"shields_crit\"/>"<<endl;
    }
}

void sediment_f::name_vtu_bedshear(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    if(p->P79==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"tau_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"tau_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==2)
    {
    result<<"<DataArray type=\"Float32\" Name=\"shearvel_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"shearvel_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==3)
    {
    result<<"<DataArray type=\"Float32\" Name=\"shields_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"shields_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
}

void sediment_f::offset_vtu_bedshear(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}

void sediment_f::print_3D_parameter1(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
    
    // alpha
    pgc->gcsl_start4(p,s->alpha,1);
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->alpha));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    
    // teta
    pgc->gcsl_start4(p,s->teta,1);
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->teta));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    // gamma
    pgc->gcsl_start4(p,s->gamma,1);
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->gamma));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    // beta
    pgc->gcsl_start4(p,s->beta,1);
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->beta));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    // phi
    pgc->gcsl_start4(p,s->phi,1);
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->phi));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
}

void sediment_f::name_pvtu_parameter1(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_alpha\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_teta\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_gamma\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_beta\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_phi\"/>"<<endl;
}

void sediment_f::name_vtu_parameter1(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_alpha\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_teta\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_gamma\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_beta\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
}

void sediment_f::offset_vtu_parameter1(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}


void sediment_f::print_3D_parameter2(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
	
    
    // dh,reduce,threshold,slideflag
    
    // dh
    pgc->gcsl_start4(p,s->vz,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->vz));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // reduce
    pgc->gcsl_start4(p,s->reduce,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->reduce));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // threshold
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=0.0;
    
    if(s->tau_eff(i,j)>s->tau_crit(i,j))
    ffn=1.0;
	result.write((char*)&ffn, sizeof (float));
	}
    
    // slideflag
    pgc->gcsl_start4(p,s->slideflag,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->slideflag));
	result.write((char*)&ffn, sizeof (float));
	}
}

void sediment_f::name_pvtu_parameter2(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_dh\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_reduce\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_threshold\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_slideflag\"/>"<<endl;
}

void sediment_f::name_vtu_parameter2(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_dh\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_reduce\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_threshold\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_slideflag\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_vtu_parameter2(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}
