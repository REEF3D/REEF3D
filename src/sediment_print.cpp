/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

void sediment_f::name_pvtu_bedload(lexer *p, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_qbe\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_qb\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_cbe\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_cb\"/>"<<endl;
}

void sediment_f::name_vtu_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_qbe\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_qb\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_cbe\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_cb\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_vtp_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_vtu_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
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

void sediment_f::print_3D_bedload(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;

	// qbe
    pgc->gcsl_start4(p,s->qbe,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->qbe));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // qb
    pgc->gcsl_start4(p,s->qb,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->qb));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // cbe
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->cbe));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // cb
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->cb));
	result.write((char*)&ffn, sizeof (float));
	}
    
}

void sediment_f::print_2D_bedload(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;

	// qbe
    pgc->gcsl_start4(p,s->qbe,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->qbe));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // qb
    pgc->gcsl_start4(p,s->qb,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->qb));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // cbe
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->cbe));
	result.write((char*)&ffn, sizeof (float));
	}
    
}

void sediment_f::name_pvtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result)
{
    if(p->P79==1)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_tau_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_tau_crit\"/>"<<endl;
    }
    
    if(p->P79==2)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shearvel_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shearvel_crit\"/>"<<endl;
    }
    
    if(p->P79==3)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shields_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shields_crit\"/>"<<endl;
    }
}

void sediment_f::name_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    if(p->P79==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==2)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==3)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
}

void sediment_f::offset_vtp_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}

void sediment_f::print_2D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
    double tau_eff,shearvel_eff,shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
    
    if(p->P79==1)
    {
        
	// tau_eff
    pgc->gcsl_start4(p,s->tau_eff,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->tau_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // tau_crit
    pgc->gcsl_start4(p,s->tau_crit,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->tau_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    
    if(p->P79==2)
    {
	// shearvel_eff
    pgc->gcsl_start4(p,s->shearvel_eff,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // shearvel_crit
    pgc->gcsl_start4(p,s->shearvel_crit,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    
    if(p->P79==3)
    {
	// shields_eff
    pgc->gcsl_start4(p,s->shields_eff,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shields_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // shields_crit
    pgc->gcsl_start4(p,s->shields_crit,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shields_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    
    }
}

void sediment_f::print_3D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
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

void sediment_f::print_2D_parameter1(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
    
    // alpha
    pgc->gcsl_start4(p,s->alpha,1);
	
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->alpha));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    
    // teta
    pgc->gcsl_start4(p,s->teta,1);
	
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->teta));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    // gamma
    pgc->gcsl_start4(p,s->gamma,1);
	
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->gamma));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    // beta
    pgc->gcsl_start4(p,s->beta,1);
	
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->beta));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
    
    // phi
    pgc->gcsl_start4(p,s->phi,1);
	
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->phi));
    ffn*=(180.0/PI);
	result.write((char*)&ffn, sizeof (float));
	}
}

void sediment_f::print_3D_parameter1(lexer* p, ghostcell *pgc, ofstream &result)
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

void sediment_f::name_pvtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_alpha\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_teta\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_gamma\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_beta\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_phi\"/>"<<endl;
}

void sediment_f::name_vtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
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

void sediment_f::offset_vtp_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_vtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
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


void sediment_f::print_2D_parameter2(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
	
    // dh,bedch,reduce,threshold,slideflag
    
    // dh
    pgc->gcsl_start4(p,s->vz,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->vz));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // bedchange
    pgc->gcsl_start4(p,s->bedch,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->bedch));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // reduce
    pgc->gcsl_start4(p,s->reduce,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->reduce));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // threshold
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=0.0;
    
    if(s->tau_eff(i,j)>s->tau_crit(i,j))
    ffn=1.0;
	result.write((char*)&ffn, sizeof (float));
	}
    
    // slideflag
    pgc->gcsl_start4(p,s->slideflag,1);
    
	iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->slideflag));
	result.write((char*)&ffn, sizeof (float));
	}
}

void sediment_f::print_3D_parameter2(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
	
    // dh,bedch,reduce,threshold,slideflag
    
    // dh
    pgc->gcsl_start4(p,s->vz,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->vz));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // bedchange
    pgc->gcsl_start4(p,s->bedch,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->bedch));
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

void sediment_f::name_pvtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_dh\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_bedch\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_reduce\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_threshold\"/>"<<endl;
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_slideflag\"/>"<<endl;
}

void sediment_f::name_vtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_dh\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_bedch\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_reduce\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_threshold\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_slideflag\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_vtp_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_vtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
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
